"""Duplicate and mass grouping utilities."""

from __future__ import annotations

import hashlib
from collections import defaultdict


def group_ids(records: list[dict]) -> None:
    formula_groups: dict[str, list[dict]] = defaultdict(list)
    exact_groups: dict[str, list[dict]] = defaultdict(list)
    nominal_groups: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        formula_groups[rec["molecular_formula"]].append(rec)
        exact_groups[f"{rec['neutral_monoisotopic_mass']:.5f}"].append(rec)
        nominal_groups[str(round(rec["neutral_monoisotopic_mass"]))].append(rec)
    for groups, field in [
        (formula_groups, "same_formula_group_id"),
        (exact_groups, "same_exact_mass_group_id"),
        (nominal_groups, "same_nominal_mass_group_id"),
    ]:
        for key, members in groups.items():
            digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:10]
            for rec in members:
                rec[field] = digest


def assign_structure_duplicate_groups(records: list[dict]) -> None:
    groups: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        groups[rec["structure_id"]].append(rec)
    for structure_id, members in groups.items():
        count = len(members)
        group_id = structure_id if count > 1 else ""
        for rec in members:
            rec["final_structure_duplicate_group_id"] = group_id
            rec["final_structure_duplicate_count"] = count
            rec["is_duplicate_final_structure"] = count > 1


def duplicate_summary(records: list[dict]) -> list[dict]:
    by_smiles: dict[str, list[dict]] = defaultdict(list)
    by_formula: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        by_smiles[rec["canonical_smiles"]].append(rec)
        by_formula[rec["molecular_formula"]].append(rec)
    rows = []
    for kind, groups in [("canonical_smiles", by_smiles), ("formula", by_formula)]:
        for key, members in groups.items():
            if len(members) > 1:
                rows.append({"group_type": kind, "group_key": key, "count": len(members), "candidate_ids": ";".join(m["candidate_id"] for m in members)})
    return rows
