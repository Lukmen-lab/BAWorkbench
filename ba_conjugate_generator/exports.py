"""Export tables, JSON, and docking snippets."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .grouping import duplicate_summary


FIELDNAMES = [
    "candidate_id", "structure_id", "candidate_name", "candidate_category", "parent_core_abbrev", "parent_core_name",
    "parent_core_smiles", "partner_name", "partner_original_label", "partner_normalized_name", "partner_smiles",
    "partner_attachment_site", "partner_attachment_atom_index",
    "reaction_channel", "conjugation_site", "modified_position", "validated_actual_modified_position",
    "site_validation_passed", "site_validation_message", "bond_type", "isomeric_smiles", "canonical_smiles",
    "core_atom_indices", "partner_atom_indices", "highlight_atom_indices", "highlight_bond_indices",
    "molecular_formula", "neutral_monoisotopic_mass", "formal_charge", "expected_charge_state_at_physiologic_pH",
    "mz_M_minus_H", "mz_M_plus_H", "mz_M_plus_Na", "mz_M_plus_NH4", "mz_M_plus_K", "mass_shift_from_parent_core",
    "same_formula_group_id", "same_exact_mass_group_id", "same_nominal_mass_group_id", "docking_ready",
    "docking_status", "docking_notes", "MS_priority", "MS_notes", "biological_plausibility", "chemical_plausibility",
    "warning_flags", "manual_review_required", "generation_status", "inclusion_reason", "exclusion_reason",
    "final_structure_duplicate_group_id", "final_structure_duplicate_count", "is_duplicate_final_structure",
    "original_candidate_id", "original_isomeric_smiles", "detected_original_site", "correction_status",
]


def _write_table(path: Path, rows: list[dict], delimiter: str = "\t", fieldnames: list[str] | None = None) -> None:
    names = fieldnames or FIELDNAMES
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=names, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def export_outputs(output: str | Path, candidates: list[dict], invalid: list[dict], excluded: list[dict]) -> None:
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    _write_table(out / "candidates_full.tsv", candidates, "\t")
    _write_table(out / "candidates_full.csv", candidates, ",")
    (out / "candidates_metadata.json").write_text(json.dumps(candidates, indent=2), encoding="utf-8")
    _write_table(out / "smiles_for_docking.tsv", candidates, "\t", ["candidate_id", "isomeric_smiles", "candidate_name"])
    manual = [row for row in candidates if str(row.get("manual_review_required", "")).strip().lower() in {"true", "1", "yes"}]
    _write_table(out / "manual_review_candidates.tsv", manual, "\t")
    _write_table(out / "invalid_candidates.tsv", invalid, "\t", sorted({k for row in invalid for k in row}) if invalid else ["partner_name", "reason"])
    _write_table(out / "excluded_partners.tsv", excluded, "\t", sorted({k for row in excluded for k in row}) if excluded else ["partner_name", "exclusion_reason"])
    _write_table(out / "duplicate_or_isomeric_groups.tsv", duplicate_summary(candidates), "\t", ["group_type", "group_key", "count", "candidate_ids"])
    summary = {
        "candidate_count": len(candidates),
        "manual_review_count": len(manual),
        "invalid_count": len(invalid),
        "excluded_partner_count": len(excluded),
        "categories": {},
        "reaction_channels": {},
    }
    for row in candidates:
        summary["categories"][row["candidate_category"]] = summary["categories"].get(row["candidate_category"], 0) + 1
        summary["reaction_channels"][row["reaction_channel"]] = summary["reaction_channels"].get(row["reaction_channel"], 0) + 1
    (out / "generation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    lig = out / "boltz_ligands"
    lig.mkdir(exist_ok=True)
    for row in candidates[:1000]:
        text = f"ligand:\n  id: {row['candidate_id']}\n  smiles: \"{row['isomeric_smiles']}\"\n  name: \"{row['candidate_name']}\"\n"
        (lig / f"{row['candidate_id']}.yaml").write_text(text, encoding="utf-8")
