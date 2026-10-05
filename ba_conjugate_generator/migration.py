"""Audit and repair positional site assignments in previous candidate exports."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from rdkit import Chem

from .cores import CORE_REGISTRY, BileAcidCore
from .exports import export_outputs
from .generator import candidate_identifiers
from .grouping import assign_structure_duplicate_groups, group_ids
from .mass import calculate_mass_properties
from .site_validation import clear_atom_maps, map_core_carbons, validate_modified_site


STEROID_O_CHANNELS = {
    "BA_O_acyl",
    "BA_O_glycoside",
    "BA_O_sulfation",
    "BA_O_phosphorylation",
    "reactive_carbonyl",
}


def _read_rows(path: str | Path) -> list[dict]:
    source = Path(path)
    delimiter = "," if source.suffix.lower() == ".csv" else "\t"
    with source.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _external_oxygen_bond(product: Chem.Mol, carbon_idx: int, core_carbons: set[int]) -> tuple[int, int, int] | None:
    carbon = product.GetAtomWithIdx(carbon_idx)
    for carbon_bond in carbon.GetBonds():
        oxygen = carbon_bond.GetOtherAtom(carbon)
        if oxygen.GetAtomicNum() != 8 or carbon_bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        external_bonds = [
            bond
            for bond in oxygen.GetBonds()
            if bond.GetOtherAtom(oxygen).GetIdx() != carbon_idx
            and bond.GetOtherAtom(oxygen).GetIdx() not in core_carbons
        ]
        if len(external_bonds) == 1:
            bond = external_bonds[0]
            return oxygen.GetIdx(), bond.GetOtherAtom(oxygen).GetIdx(), bond.GetIdx()
    return None


def _target_hydroxyl_oxygen(product: Chem.Mol, carbon_idx: int) -> int:
    carbon = product.GetAtomWithIdx(carbon_idx)
    oxygens = [
        bond.GetOtherAtom(carbon).GetIdx()
        for bond in carbon.GetBonds()
        if bond.GetBondType() == Chem.BondType.SINGLE and bond.GetOtherAtom(carbon).GetAtomicNum() == 8
    ]
    if len(oxygens) != 1:
        raise ValueError(f"Expected one steroid oxygen on mapped carbon atom {carbon_idx}; found {len(oxygens)}")
    return oxygens[0]


def repair_product_site(core: BileAcidCore, product: Chem.Mol, actual_position: str, requested_position: str) -> Chem.Mol:
    requested_site = core.site(requested_position)
    if requested_site is None or requested_site.site_kind != "steroid_hydroxyl":
        raise ValueError(f"{core.abbrev} has no steroid hydroxyl at {requested_position}")
    carbon_map, error = map_core_carbons(core, product)
    if error:
        raise ValueError(error)
    core_carbons = set(carbon_map.values())
    actual_carbon = carbon_map[24] if actual_position == "C24" else carbon_map[int(actual_position[1:])]
    attachment = _external_oxygen_bond(product, actual_carbon, core_carbons)
    if attachment is None:
        raise ValueError(f"Could not locate external substituent bond at detected position {actual_position}")
    source_oxygen, partner_atom, _old_bond = attachment
    target_oxygen = _target_hydroxyl_oxygen(product, carbon_map[requested_site.carbon_map_num])
    if source_oxygen == target_oxygen:
        return Chem.Mol(product)
    rw = Chem.RWMol(product)
    rw.RemoveBond(source_oxygen, partner_atom)
    source_atom = rw.GetAtomWithIdx(source_oxygen)
    source_atom.SetNumExplicitHs(0)
    source_atom.SetNoImplicit(False)
    target_atom = rw.GetAtomWithIdx(target_oxygen)
    target_atom.SetNumExplicitHs(0)
    target_atom.SetNoImplicit(True)
    rw.AddBond(target_oxygen, partner_atom, Chem.BondType.SINGLE)
    repaired = rw.GetMol()
    Chem.SanitizeMol(repaired)
    Chem.AssignStereochemistry(repaired, cleanIt=True, force=True)
    return repaired


def _apply_corrected_structure(row: dict, mol: Chem.Mol, validation) -> dict:
    clean = clear_atom_maps(mol)
    smiles = Chem.MolToSmiles(clean, isomericSmiles=True, canonical=True)
    canonical = Chem.MolToSmiles(clean, isomericSmiles=False, canonical=True)
    props = calculate_mass_properties(clean)
    candidate_id, structure_id = candidate_identifiers(
        smiles,
        row["parent_core_abbrev"],
        row.get("partner_normalized_name", ""),
        row["reaction_channel"],
        row["modified_position"],
    )
    corrected = dict(row)
    corrected.update(
        {
            "candidate_id": candidate_id,
            "structure_id": structure_id,
            "isomeric_smiles": smiles,
            "canonical_smiles": canonical,
            "molecular_formula": props.formula,
            "neutral_monoisotopic_mass": round(props.exact_mass, 6),
            "formal_charge": props.formal_charge,
            "mz_M_minus_H": round(props.mz_m_minus_h, 6),
            "mz_M_plus_H": round(props.mz_m_plus_h, 6),
            "mz_M_plus_Na": round(props.mz_m_plus_na, 6),
            "mz_M_plus_NH4": round(props.mz_m_plus_nh4, 6),
            "mz_M_plus_K": round(props.mz_m_plus_k, 6),
            "validated_actual_modified_position": validation.actual_position,
            "site_validation_passed": validation.passed,
            "site_validation_message": validation.message,
            "core_atom_indices": ";".join(map(str, validation.core_atom_indices)),
            "partner_atom_indices": ";".join(map(str, validation.partner_atom_indices)),
            "highlight_atom_indices": ";".join(map(str, validation.modified_atom_indices)),
            "highlight_bond_indices": ";".join(map(str, validation.new_bond_indices)),
        }
    )
    if not validation.passed:
        corrected["docking_ready"] = False
        corrected["docking_status"] = "site_validation_failed"
    return corrected


def migrate_previous_output(input_path: str | Path, output_dir: str | Path) -> dict[str, object]:
    source_rows = _read_rows(input_path)
    corrected_rows: list[dict] = []
    audit_rows: list[dict] = []
    mismatch_count = 0
    corrected_count = 0
    failed_count = 0
    for row in source_rows:
        original_id = row["candidate_id"]
        original_smiles = row["isomeric_smiles"]
        core = CORE_REGISTRY[row["parent_core_abbrev"]]
        original_mol = Chem.MolFromSmiles(original_smiles)
        if original_mol is None:
            original_validation = None
            status = "failed_invalid_original_smiles"
            message = "Original isomeric SMILES could not be parsed."
            corrected = dict(row)
            failed_count += 1
        else:
            original_validation = validate_modified_site(core, original_mol, row["modified_position"])
            if not original_validation.passed:
                mismatch_count += 1
            try:
                if original_validation.passed:
                    repaired_mol = original_mol
                    status = "unchanged_valid"
                elif row["reaction_channel"] in STEROID_O_CHANNELS:
                    repaired_mol = repair_product_site(core, original_mol, original_validation.actual_position, row["modified_position"])
                    status = "corrected"
                else:
                    raise ValueError("Mismatch is not in a supported steroid-O migration channel")
                final_validation = validate_modified_site(core, repaired_mol, row["modified_position"])
                corrected = _apply_corrected_structure(row, repaired_mol, final_validation)
                if status == "corrected" and final_validation.passed:
                    corrected_count += 1
                if not final_validation.passed:
                    status = "failed_post_validation"
                    failed_count += 1
                message = final_validation.message
            except Exception as exc:
                corrected = dict(row)
                corrected["docking_ready"] = False
                corrected["docking_status"] = "site_validation_failed"
                status = "failed_repair"
                message = str(exc)
                failed_count += 1
        corrected["original_candidate_id"] = original_id
        corrected["original_isomeric_smiles"] = original_smiles
        corrected["detected_original_site"] = original_validation.actual_position if original_validation else "unresolved"
        corrected["correction_status"] = status
        corrected_rows.append(corrected)
        old_mass = float(row["neutral_monoisotopic_mass"])
        new_mass = float(corrected.get("neutral_monoisotopic_mass", old_mass))
        audit_rows.append(
            {
                "original_candidate_id": original_id,
                "original_isomeric_smiles": original_smiles,
                "requested_modified_position": row.get("modified_position", ""),
                "detected_original_site": corrected["detected_original_site"],
                "corrected_candidate_id": corrected.get("candidate_id", ""),
                "corrected_isomeric_smiles": corrected.get("isomeric_smiles", ""),
                "correction_status": status,
                "site_validation_message": message,
                "formula_unchanged": str(row.get("molecular_formula", "")) == str(corrected.get("molecular_formula", "")),
                "exact_mass_unchanged": abs(old_mass - new_mass) < 1e-6,
                "formal_charge_unchanged": int(float(row.get("formal_charge", 0))) == int(float(corrected.get("formal_charge", 0))),
            }
        )
    group_ids(corrected_rows)
    assign_structure_duplicate_groups(corrected_rows)
    export_outputs(output_dir, corrected_rows, [], [])
    out = Path(output_dir)
    audit_fields = list(audit_rows[0]) if audit_rows else []
    with (out / "site_migration_audit.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)
    summary = {
        "input": str(input_path),
        "total_rows": len(source_rows),
        "mismatches_before_repair": mismatch_count,
        "corrected_rows": corrected_count,
        "failed_rows": failed_count,
        "mismatches_after_repair": sum(not _as_bool(row.get("site_validation_passed")) for row in corrected_rows),
        "formula_changes": sum(not row["formula_unchanged"] for row in audit_rows),
        "exact_mass_changes": sum(not row["exact_mass_unchanged"] for row in audit_rows),
        "formal_charge_changes": sum(not row["formal_charge_unchanged"] for row in audit_rows),
    }
    (out / "site_migration_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
