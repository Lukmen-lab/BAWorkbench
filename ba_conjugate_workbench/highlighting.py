"""Highlight heuristics for generated candidate molecules."""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

from ba_conjugate_generator.cores import CORE_REGISTRY
from ba_conjugate_generator.site_validation import map_core_carbons, validate_modified_site

from .utils import candidate_smiles


@dataclass
class HighlightSpec:
    core_atoms: list[int] = field(default_factory=list)
    partner_atoms: list[int] = field(default_factory=list)
    site_atoms: list[int] = field(default_factory=list)
    new_bonds: list[int] = field(default_factory=list)
    labels: dict[int, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def _parse_indices(value: object) -> list[int]:
    if value is None:
        return []
    out = []
    for part in str(value).replace(",", ";").split(";"):
        part = part.strip()
        if part.isdigit():
            out.append(int(part))
    return out


def _bond_between(mol: Chem.Mol, a: int, b: int) -> int | None:
    bond = mol.GetBondBetweenAtoms(a, b)
    return bond.GetIdx() if bond is not None else None


def conventional_ba_labels(row: dict) -> dict[int, str]:
    """Return standard bile-acid labels mapped to product atom indices.

    The labels are conventional BA positions, not RDKit atom numbers.
    """
    core_key = str(row.get("parent_core_abbrev", ""))
    core = CORE_REGISTRY.get(core_key)
    if core is None:
        return {}
    smiles = candidate_smiles(row)
    product = Chem.MolFromSmiles(smiles) if smiles else None
    if product is None:
        return {}
    carbon_map, error = map_core_carbons(core, product)
    if error:
        return {}
    labels: dict[int, str] = {}
    for site in core.hydroxyl_site_objects:
        carbon_idx = carbon_map[site.carbon_map_num]
        labels[carbon_idx] = site.position
        single_oxygen_neighbors = [
            bond.GetOtherAtom(product.GetAtomWithIdx(carbon_idx)).GetIdx()
            for bond in product.GetAtomWithIdx(carbon_idx).GetBonds()
            if bond.GetBondType() == Chem.BondType.SINGLE and bond.GetOtherAtom(product.GetAtomWithIdx(carbon_idx)).GetAtomicNum() == 8
        ]
        if len(single_oxygen_neighbors) == 1:
            labels[single_oxygen_neighbors[0]] = f"{site.position}-O"
    labels[carbon_map[24]] = "C24"
    return labels


def infer_highlights(row: dict) -> HighlightSpec:
    smiles = candidate_smiles(row)
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    spec = HighlightSpec()
    if mol is None:
        spec.warnings.append("Highlighting unavailable: invalid or missing SMILES.")
        return spec
    spec.core_atoms = _parse_indices(row.get("core_atom_indices"))
    spec.partner_atoms = _parse_indices(row.get("partner_atom_indices"))
    spec.site_atoms = _parse_indices(row.get("highlight_atom_indices"))
    spec.new_bonds = _parse_indices(row.get("highlight_bond_indices") or row.get("new_bond_indices"))
    core = CORE_REGISTRY.get(str(row.get("parent_core_abbrev", "")))
    if core is not None:
        validation = validate_modified_site(core, mol, str(row.get("modified_position", "")))
        if not spec.core_atoms:
            spec.core_atoms = list(validation.core_atom_indices)
        if not spec.partner_atoms:
            spec.partner_atoms = list(validation.partner_atom_indices)
        if not spec.site_atoms:
            spec.site_atoms = list(validation.modified_atom_indices)
        if not spec.new_bonds:
            spec.new_bonds = list(validation.new_bond_indices)
    if not spec.partner_atoms and spec.core_atoms:
        core = set(spec.core_atoms)
        spec.partner_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in core]
    spec.labels.update(conventional_ba_labels(row))
    if not spec.core_atoms:
        spec.warnings.append("BA core heuristic did not match; core highlighting may be incomplete.")
    return spec


def large_flexible_warning(row: dict) -> str | None:
    flags = str(row.get("warning_flags", ""))
    mass = row.get("neutral_monoisotopic_mass")
    try:
        large_mass = float(mass) > 900
    except Exception:
        large_mass = False
    if "large_flexible_warning" in flags or "polymeric_fragment_model" in flags or large_mass:
        return "This molecule is large/flexible. A single generated conformer may not represent the relevant docking conformation."
    return None
