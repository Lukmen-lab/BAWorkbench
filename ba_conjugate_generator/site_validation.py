"""Graph-based validation of the actual bile-acid modification site."""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

from .cores import BileAcidCore


@dataclass(frozen=True)
class SiteValidationResult:
    actual_position: str
    passed: bool
    message: str
    core_atom_indices: tuple[int, ...] = field(default_factory=tuple)
    partner_atom_indices: tuple[int, ...] = field(default_factory=tuple)
    modified_atom_indices: tuple[int, ...] = field(default_factory=tuple)
    new_bond_indices: tuple[int, ...] = field(default_factory=tuple)


def clear_atom_maps(mol: Chem.Mol) -> Chem.Mol:
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(clean)
    Chem.AssignStereochemistry(clean, cleanIt=True, force=True)
    return clean


def _carbon_scaffold_query(core: BileAcidCore) -> Chem.Mol:
    """Build a carbon-only cholane query carrying the stored conventional maps."""
    mapped = core.mapped_mol
    rw = Chem.RWMol()
    source_to_query: dict[int, int] = {}
    for atom in mapped.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        query_atom = Chem.AtomFromSmarts("[#6]")
        query_atom.SetAtomMapNum(atom.GetAtomMapNum())
        source_to_query[atom.GetIdx()] = rw.AddAtom(query_atom)
    for bond in mapped.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in source_to_query and end in source_to_query:
            rw.AddBond(source_to_query[begin], source_to_query[end], Chem.BondType.SINGLE)
    return rw.GetMol()


def map_core_carbons(core: BileAcidCore, product: Chem.Mol) -> tuple[dict[int, int], str]:
    query = _carbon_scaffold_query(core)
    matches = product.GetSubstructMatches(query, uniquify=True, useChirality=False)
    if len(matches) != 1:
        return {}, f"Expected one mapped {core.abbrev} carbon scaffold in product; found {len(matches)}."
    conventional_to_product: dict[int, int] = {}
    for query_idx, product_idx in enumerate(matches[0]):
        map_num = query.GetAtomWithIdx(query_idx).GetAtomMapNum()
        conventional_to_product[map_num] = product_idx
    return conventional_to_product, ""


def _external_bond_from_steroid_oxygen(
    product: Chem.Mol,
    carbon_idx: int,
    core_carbon_indices: set[int],
) -> tuple[int, int] | None:
    carbon = product.GetAtomWithIdx(carbon_idx)
    for bond in carbon.GetBonds():
        oxygen = bond.GetOtherAtom(carbon)
        if oxygen.GetAtomicNum() != 8 or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        for oxygen_bond in oxygen.GetBonds():
            external = oxygen_bond.GetOtherAtom(oxygen)
            if external.GetIdx() == carbon_idx:
                continue
            if external.GetIdx() not in core_carbon_indices:
                return oxygen.GetIdx(), oxygen_bond.GetIdx()
    return None


def _external_bond_at_c24(
    product: Chem.Mol,
    c24_idx: int,
    c23_idx: int,
    core_carbon_indices: set[int],
) -> tuple[int, int] | None:
    c24 = product.GetAtomWithIdx(c24_idx)
    for bond in c24.GetBonds():
        neighbor = bond.GetOtherAtom(c24)
        if neighbor.GetIdx() == c23_idx or bond.GetBondType() == Chem.BondType.DOUBLE:
            continue
        if neighbor.GetAtomicNum() == 7:
            return neighbor.GetIdx(), bond.GetIdx()
        if neighbor.GetAtomicNum() == 8:
            for oxygen_bond in neighbor.GetBonds():
                external = oxygen_bond.GetOtherAtom(neighbor)
                if external.GetIdx() == c24_idx:
                    continue
                if external.GetIdx() not in core_carbon_indices:
                    return neighbor.GetIdx(), oxygen_bond.GetIdx()
    return None


def validate_modified_site(core: BileAcidCore, product: Chem.Mol, requested_position: str) -> SiteValidationResult:
    carbon_map, mapping_error = map_core_carbons(core, product)
    if mapping_error:
        return SiteValidationResult("unresolved", False, mapping_error)
    core_carbons = set(carbon_map.values())
    core_atoms = set(core_carbons)
    for carbon_idx in core_carbons:
        carbon = product.GetAtomWithIdx(carbon_idx)
        core_atoms.update(neighbor.GetIdx() for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 8)
    detected: list[tuple[str, int, int]] = []
    for site in core.hydroxyl_site_objects:
        external = _external_bond_from_steroid_oxygen(product, carbon_map[site.carbon_map_num], core_carbons)
        if external is not None:
            detected.append((site.position, external[0], external[1]))
    c24_external = _external_bond_at_c24(product, carbon_map[24], carbon_map[23], core_carbons)
    if c24_external is not None:
        detected.append(("C24", c24_external[0], c24_external[1]))
    if len(detected) != 1:
        actual = "none" if not detected else ";".join(position for position, _, _ in detected)
        return SiteValidationResult(
            actual,
            False,
            f"Expected exactly one modified BA position; detected {actual}.",
            tuple(sorted(core_atoms)),
            tuple(sorted(set(range(product.GetNumAtoms())) - core_atoms)),
        )
    actual, modified_atom, new_bond = detected[0]
    passed = actual == requested_position
    message = (
        f"Validated requested BA position {requested_position}."
        if passed
        else f"Requested BA position {requested_position}, but graph validation detected {actual}."
    )
    partner_atoms = tuple(sorted(set(range(product.GetNumAtoms())) - core_atoms))
    return SiteValidationResult(
        actual,
        passed,
        message,
        tuple(sorted(core_atoms)),
        partner_atoms,
        (carbon_map[24] if actual == "C24" else carbon_map[int(actual[1:])], modified_atom),
        (new_bond,),
    )
