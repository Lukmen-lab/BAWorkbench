"""Table presentation helpers."""

from __future__ import annotations

import pandas as pd


DEFAULT_COLUMNS = [
    "candidate_id",
    "candidate_name",
    "parent_core_abbrev",
    "candidate_category",
    "partner_name",
    "partner_attachment_site",
    "reaction_channel",
    "conjugation_site",
    "modified_position",
    "validated_actual_modified_position",
    "site_validation_passed",
    "molecular_formula",
    "neutral_monoisotopic_mass",
    "mz_M_minus_H",
    "mz_M_plus_H",
    "docking_status",
    "warning_flags",
    "curation_status",
]

OPTIONAL_COLUMNS = [
    "isomeric_smiles",
    "canonical_smiles",
    "formal_charge",
    "mass_shift_from_parent_core",
    "same_formula_group_id",
    "same_exact_mass_group_id",
    "same_nominal_mass_group_id",
    "biological_plausibility",
    "chemical_plausibility",
    "manual_review_required",
    "MS_notes",
    "docking_notes",
]

DETAIL_COLUMNS = [
    "candidate_id",
    "candidate_name",
    "parent_core_abbrev",
    "parent_core_name",
    "partner_name",
    "partner_original_label",
    "partner_normalized_name",
    "partner_attachment_site",
    "partner_attachment_atom_index",
    "reaction_channel",
    "conjugation_site",
    "modified_position",
    "validated_actual_modified_position",
    "site_validation_passed",
    "site_validation_message",
    "bond_type",
    "isomeric_smiles",
    "canonical_smiles",
    "molecular_formula",
    "neutral_monoisotopic_mass",
    "formal_charge",
    "expected_charge_state_at_physiologic_pH",
    "mz_M_minus_H",
    "mz_M_plus_H",
    "mz_M_plus_Na",
    "mz_M_plus_NH4",
    "mz_M_plus_K",
    "mass_shift_from_parent_core",
    "docking_status",
    "docking_notes",
    "MS_priority",
    "MS_notes",
    "warning_flags",
    "manual_review_required",
    "biological_plausibility",
    "chemical_plausibility",
    "curation_status",
    "curation_notes",
]


def display_columns(df: pd.DataFrame, include_optional: list[str] | None = None) -> list[str]:
    include_optional = include_optional or []
    cols = [c for c in DEFAULT_COLUMNS if c in df.columns]
    cols.extend(c for c in include_optional if c in df.columns and c not in cols)
    return cols


def candidate_id_from_selection(table: pd.DataFrame, selected_rows: list[int]) -> str | None:
    """Resolve a native dataframe row-selection event to a candidate ID."""
    if not selected_rows or "candidate_id" not in table.columns:
        return None
    row_index = selected_rows[-1]
    if row_index < 0 or row_index >= len(table):
        return None
    return str(table.iloc[row_index]["candidate_id"])
