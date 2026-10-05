"""Candidate generation orchestration."""

from __future__ import annotations

import hashlib
from dataclasses import asdict

from rdkit import Chem

from .config import GeneratorConfig
from .cores import BileAcidCore
from .grouping import assign_structure_duplicate_groups, group_ids
from .mass import calculate_mass_properties
from .partners import Partner
from .reactions import (
    ReactionProduct,
    ba_o_acyl,
    c24_amides,
    c24_esters,
    glucan_fragment_products,
    glycosides,
    reactive_carbonyl_adducts,
    sulfate_phosphate_esters,
)
from .site_validation import clear_atom_maps, validate_modified_site


PRIMARY_CHANNELS = {"C24_amide", "C24_ester", "BA_O_acyl", "BA_O_glycoside"}


def candidate_identifiers(smiles: str, core_abbrev: str, partner_normalized_name: str, reaction_channel: str, position: str) -> tuple[str, str]:
    structure_id = f"BASTR_{hashlib.sha256(smiles.encode('utf-8')).hexdigest()[:24]}"
    token = hashlib.sha1(f"{structure_id}|{core_abbrev}|{partner_normalized_name}|{reaction_channel}|{position}".encode("utf-8")).hexdigest()[:12]
    return f"{core_abbrev}_{token}", structure_id


def _category_for_partner(base: str, partner: Partner) -> str:
    if partner.source_category == "vitamin" and base not in {"BA_reactive_carbonyl_adduct"}:
        return "BA_vitamin_conjugate"
    if partner.source_category == "nucleobase":
        return "BA_nucleobase_conjugate"
    if partner.source_category == "indole_metabolite":
        return "BA_indole_metabolite_conjugate"
    return base


def unmodified_record(core: BileAcidCore, config: GeneratorConfig) -> dict:
    """Build a candidate row for an unmodified parent bile acid."""
    mol = core.mol
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    can = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    mass = calculate_mass_properties(mol)
    candidate_id, structure_id = candidate_identifiers(smiles, core.abbrev, "none", "unmodified", "none")
    # Parent BA references are validated inputs, not exploratory generated
    # products. Keep them available even when manual-review candidates are
    # filtered from a manual-review generation run.
    manual = False
    return {
        "candidate_id": candidate_id,
        "structure_id": structure_id,
        "candidate_name": f"{core.name} (unmodified)",
        "candidate_category": "unmodified_BA",
        "parent_core_abbrev": core.abbrev,
        "parent_core_name": core.name,
        "parent_core_smiles": core.smiles,
        "partner_name": "none",
        "partner_original_label": "none",
        "partner_normalized_name": "none",
        "partner_smiles": "",
        "partner_attachment_site": "none",
        "partner_attachment_atom_index": "",
        "reaction_channel": "unmodified",
        "conjugation_site": "none",
        "modified_position": "none",
        "validated_actual_modified_position": "none",
        "site_validation_passed": True,
        "site_validation_message": "Unmodified parent BA; no modification site is expected.",
        "core_atom_indices": ";".join(str(atom.GetIdx()) for atom in mol.GetAtoms()),
        "partner_atom_indices": "",
        "highlight_atom_indices": "",
        "highlight_bond_indices": "",
        "bond_type": "none",
        "isomeric_smiles": smiles,
        "canonical_smiles": can,
        "molecular_formula": mass.formula,
        "neutral_monoisotopic_mass": round(mass.exact_mass, 6),
        "formal_charge": mass.formal_charge,
        "expected_charge_state_at_physiologic_pH": "neutral_or_anionic",
        "mz_M_minus_H": round(mass.mz_m_minus_h, 6),
        "mz_M_plus_H": round(mass.mz_m_plus_h, 6),
        "mz_M_plus_Na": round(mass.mz_m_plus_na, 6),
        "mz_M_plus_NH4": round(mass.mz_m_plus_nh4, 6),
        "mz_M_plus_K": round(mass.mz_m_plus_k, 6),
        "mass_shift_from_parent_core": 0.0,
        "same_formula_group_id": "",
        "same_exact_mass_group_id": "",
        "same_nominal_mass_group_id": "",
        "docking_ready": not manual,
        "docking_status": "manual_review_required" if manual else "ready",
        "docking_notes": "Unmodified parent BA included as a reference ligand.",
        "MS_priority": "high",
        "MS_notes": "Use as the unmodified parent BA reference for exact-mass and LC-MS/MS comparisons.",
        "biological_plausibility": "reference_parent",
        "chemical_plausibility": "known_parent_structure",
        "warning_flags": "",
        "manual_review_required": manual,
        "generation_status": "included_unmodified",
        "inclusion_reason": "Selected unmodified parent BA included as a reference candidate.",
        "exclusion_reason": "",
    }


def _record(core: BileAcidCore, partner: Partner, product: ReactionProduct, config: GeneratorConfig) -> dict:
    Chem.SanitizeMol(product.mol)
    Chem.AssignStereochemistry(product.mol, cleanIt=True, force=True)
    validation = validate_modified_site(core, product.mol, product.site.position)
    clean_mol = clear_atom_maps(product.mol)
    smiles = Chem.MolToSmiles(clean_mol, isomericSmiles=True, canonical=True)
    can = Chem.MolToSmiles(clean_mol, isomericSmiles=False, canonical=True)
    mass = calculate_mass_properties(clean_mol)
    parent_mass = calculate_mass_properties(core.mol).exact_mass
    candidate_id, structure_id = candidate_identifiers(smiles, core.abbrev, partner.normalized_name, product.reaction_channel, product.site.position)
    warnings = sorted(set(product.warning_flags + partner.warning_flags))
    manual = "manual_review_required" in warnings or config.mode == "manual_review"
    docking_ready = validation.passed and not manual and "large_flexible_warning" not in warnings
    status = "site_validation_failed" if not validation.passed else ("manual_review_required" if manual else ("large_flexible_review" if "large_flexible_warning" in warnings else "ready"))
    category = _category_for_partner(product.candidate_category, partner)
    attachment_name = f" ({product.partner_attachment_site})" if product.partner_attachment_site else ""
    return {
        "candidate_id": candidate_id,
        "structure_id": structure_id,
        "candidate_name": f"{core.abbrev} {product.modified_position} {partner.partner_name} {product.bond_type}{attachment_name}",
        "candidate_category": category,
        "parent_core_abbrev": core.abbrev,
        "parent_core_name": core.name,
        "parent_core_smiles": core.smiles,
        "partner_name": partner.partner_name,
        "partner_original_label": ";".join(partner.original_labels or [partner.partner_name]),
        "partner_normalized_name": partner.normalized_name,
        "partner_smiles": partner.smiles or "",
        "partner_attachment_site": product.partner_attachment_site,
        "partner_attachment_atom_index": (
            product.partner_attachment_atom_index if product.partner_attachment_atom_index is not None else ""
        ),
        "reaction_channel": product.reaction_channel,
        "conjugation_site": product.conjugation_site,
        "modified_position": product.modified_position,
        "validated_actual_modified_position": validation.actual_position,
        "site_validation_passed": validation.passed,
        "site_validation_message": validation.message,
        "core_atom_indices": ";".join(map(str, validation.core_atom_indices)),
        "partner_atom_indices": ";".join(map(str, validation.partner_atom_indices)),
        "highlight_atom_indices": ";".join(map(str, validation.modified_atom_indices)),
        "highlight_bond_indices": ";".join(map(str, validation.new_bond_indices)),
        "bond_type": product.bond_type,
        "isomeric_smiles": smiles,
        "canonical_smiles": can,
        "molecular_formula": mass.formula,
        "neutral_monoisotopic_mass": round(mass.exact_mass, 6),
        "formal_charge": mass.formal_charge,
        "expected_charge_state_at_physiologic_pH": "neutral_or_anionic",
        "mz_M_minus_H": round(mass.mz_m_minus_h, 6),
        "mz_M_plus_H": round(mass.mz_m_plus_h, 6),
        "mz_M_plus_Na": round(mass.mz_m_plus_na, 6),
        "mz_M_plus_NH4": round(mass.mz_m_plus_nh4, 6),
        "mz_M_plus_K": round(mass.mz_m_plus_k, 6),
        "mass_shift_from_parent_core": round(mass.exact_mass - parent_mass, 6),
        "same_formula_group_id": "",
        "same_exact_mass_group_id": "",
        "same_nominal_mass_group_id": "",
        "docking_ready": docking_ready,
        "docking_status": status,
        "docking_notes": product.generation_note,
        "MS_priority": "high" if product.reaction_channel in {"C24_amide", "C24_ester"} else ("manual_review" if manual else "medium"),
        "MS_notes": "Use exact mass, formula, expected adducts, and warning flags for LC-MS/MS follow-up.",
        "biological_plausibility": partner.plausibility_class,
        "chemical_plausibility": "exploratory" if manual else "chemically_plausible",
        "warning_flags": ";".join(warnings),
        "manual_review_required": manual,
        "generation_status": "generated",
        "inclusion_reason": "Chemically defined covalent BA-derived candidate generated by RDKit-edit reaction channel.",
        "exclusion_reason": "",
    }


def generate_candidates(cores: list[BileAcidCore], partners: list[Partner], config: GeneratorConfig) -> tuple[list[dict], list[dict], list[dict]]:
    candidates: list[dict] = []
    invalid: list[dict] = []
    excluded: list[dict] = []
    for core in cores:
        try:
            candidates.append(unmodified_record(core, config))
        except Exception as exc:
            invalid.append({"partner_name": "none", "core": core.abbrev, "reason": f"Unmodified BA: {exc}"})
    for partner in partners:
        if partner.is_excluded_context_only:
            excluded.append({**asdict(partner), "exclusion_reason": "Context-only salt/inorganic/metal-complex/undefined input; no covalent BA candidate generated."})
            continue
        if partner.smiles and partner.mol is None:
            invalid.append({"partner_name": partner.partner_name, "partner_smiles": partner.smiles, "reason": "RDKit failed to parse or sanitize partner SMILES"})
            continue
        if not partner.allowed_reaction_channels:
            excluded.append({**asdict(partner), "exclusion_reason": "No chemically meaningful covalent BA reaction channel in prototype."})
            continue
        for core in cores:
            products: list[ReactionProduct] = []
            products.extend(c24_amides(core, partner))
            products.extend(c24_esters(core, partner))
            products.extend(ba_o_acyl(core, partner))
            if config.mode != "primary_docking" or partner.source_category == "saccharide":
                products.extend(glycosides(core, partner))
            if config.mode in {"broad_covalent", "manual_review"}:
                products.extend(reactive_carbonyl_adducts(core, partner))
            if config.mode == "primary_docking":
                products = [p for p in products if p.reaction_channel in PRIMARY_CHANNELS and p.candidate_category != "BA_polymeric_fragment_model"]
            for product in products:
                try:
                    candidates.append(_record(core, partner, product, config))
                except Exception as exc:
                    invalid.append({"partner_name": partner.partner_name, "core": core.abbrev, "reason": str(exc)})
        if "glucan_fragments" in partner.allowed_reaction_channels and config.mode in {"broad_covalent", "manual_review"}:
            for core in cores:
                for product in glucan_fragment_products(core, config.max_saccharide_dp):
                    candidates.append(_record(core, partner, product, config))
    if config.mode in {"broad_covalent", "manual_review"}:
        pseudo = Partner("inorganic sulfate/phosphate donor", "covalent_sulfate_phosphate", smiles="", allowed_reaction_channels=[])
        for core in cores:
            for product in sulfate_phosphate_esters(core):
                candidates.append(_record(core, pseudo, product, config))
    unique_candidates: list[dict] = []
    candidates_by_id: dict[str, dict] = {}
    for candidate in candidates:
        candidate_id = candidate["candidate_id"]
        existing = candidates_by_id.get(candidate_id)
        if existing is None:
            candidates_by_id[candidate_id] = candidate
            unique_candidates.append(candidate)
            continue
        merged_warnings = sorted(
            {
                warning
                for value in (existing.get("warning_flags", ""), candidate.get("warning_flags", ""))
                for warning in str(value).split(";")
                if warning
            }
        )
        existing["warning_flags"] = ";".join(merged_warnings)
    candidates = unique_candidates
    group_ids(candidates)
    assign_structure_duplicate_groups(candidates)
    return candidates, invalid, excluded
