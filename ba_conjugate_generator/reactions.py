"""Covalent reaction channels implemented with RDKit molecule edits."""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

from .cores import BileAcidCore, BileAcidSite
from .partners import Partner
from .validation import canonical_smiles


CARBOXYL = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
AMINE = Chem.MolFromSmarts("[NX3;H1,H2;!$(NC=O)]")
ALCOHOL = Chem.MolFromSmarts("[OX2H1][#6]")
CARBONYL = Chem.MolFromSmarts("[CX3]=[OX1]")


@dataclass
class ReactionProduct:
    mol: Chem.Mol
    reaction_channel: str
    candidate_category: str
    site: BileAcidSite
    bond_type: str
    warning_flags: list[str] = field(default_factory=list)
    generation_note: str = ""
    partner_attachment_site: str = ""
    partner_attachment_atom_index: int | None = None

    @property
    def conjugation_site(self) -> str:
        return self.site.conjugation_site

    @property
    def modified_position(self) -> str:
        return self.site.position


@dataclass(frozen=True)
class SiteResolution:
    sites: tuple[BileAcidSite, ...]
    skipped_reason: str = ""


def resolve_steroid_sites(core: BileAcidCore, requested_position: str | None = None) -> SiteResolution:
    if requested_position is None:
        return SiteResolution(core.hydroxyl_site_objects)
    site = core.site(requested_position)
    if site is None or site.site_kind != "steroid_hydroxyl":
        return SiteResolution((), f"{core.abbrev} has no steroid hydroxyl at {requested_position}; no substitute site was used.")
    return SiteResolution((site,))


def _mapped_atom_idx(mol: Chem.Mol, atom_map_num: int) -> int:
    indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == atom_map_num]
    if len(indices) != 1:
        raise ValueError(f"Expected one atom with map {atom_map_num}, found {len(indices)}")
    return indices[0]


def _sanitize_and_assign(mol: Chem.Mol) -> Chem.Mol:
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def _core_prepared_for_oxygen_bond(core: BileAcidCore, site: BileAcidSite) -> tuple[Chem.Mol, int]:
    mapped_core = core.mapped_mol
    oxygen_idx = _mapped_atom_idx(mapped_core, site.oxygen_map_num)
    rw = Chem.RWMol(mapped_core)
    oxygen = rw.GetAtomWithIdx(oxygen_idx)
    oxygen.SetNumExplicitHs(0)
    oxygen.SetNoImplicit(True)
    return rw.GetMol(), oxygen_idx


def _remove_atom(mol: Chem.Mol, atom_idx: int) -> tuple[Chem.Mol, dict[int, int]]:
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(atom_idx)
    mapping = {old: old - (1 if old > atom_idx else 0) for old in range(mol.GetNumAtoms()) if old != atom_idx}
    out = rw.GetMol()
    return _sanitize_and_assign(out), mapping


def _combine_and_bond(left: Chem.Mol, left_idx: int, right: Chem.Mol, right_idx: int, order=Chem.BondType.SINGLE) -> Chem.Mol | None:
    combo = Chem.CombineMols(left, right)
    rw = Chem.RWMol(combo)
    rw.AddBond(left_idx, left.GetNumAtoms() + right_idx, order)
    mol = rw.GetMol()
    try:
        _sanitize_and_assign(mol)
    except Exception:
        return None
    return mol


def _carboxyl_sites(mol: Chem.Mol) -> list[tuple[int, int]]:
    sites: list[tuple[int, int]] = []
    for match in mol.GetSubstructMatches(CARBOXYL):
        carbon, _oxo, hydroxyl = match
        sites.append((carbon, hydroxyl))
    return sites


def _amine_sites(mol: Chem.Mol) -> list[int]:
    return sorted({m[0] for m in mol.GetSubstructMatches(AMINE)})


def _arginine_amide_sites(mol: Chem.Mol) -> list[tuple[int, str]]:
    """Return the alpha N and both terminal guanidino Ns.

    The chain-connected guanidino N is deliberately excluded. The two
    terminal nitrogens have different bond orders only because the neutral
    guanidine input requires one explicit tautomer; both are retained as
    separate attachment models.
    """
    guanidino_carbons = []
    for atom in mol.GetAtoms():
        nitrogen_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 7]
        if atom.GetAtomicNum() == 6 and len(nitrogen_neighbors) == 3:
            guanidino_carbons.append(atom)
    if len(guanidino_carbons) != 1:
        raise ValueError(f"Expected one arginine guanidino carbon; found {len(guanidino_carbons)}.")

    guanidino_carbon = guanidino_carbons[0]
    guanidino_nitrogens = list(guanidino_carbon.GetNeighbors())
    terminal_nitrogens = [
        nitrogen
        for nitrogen in guanidino_nitrogens
        if nitrogen.GetAtomicNum() == 7
        and len([neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetAtomicNum() > 1]) == 1
    ]
    alpha_nitrogens = [
        mol.GetAtomWithIdx(idx)
        for idx in _amine_sites(mol)
        if idx not in {nitrogen.GetIdx() for nitrogen in guanidino_nitrogens}
    ]
    if len(alpha_nitrogens) != 1 or len(terminal_nitrogens) != 2:
        raise ValueError(
            "Expected one alpha-amino N and two terminal guanidino Ns for arginine; "
            f"found {len(alpha_nitrogens)} and {len(terminal_nitrogens)}."
        )

    terminal_nitrogens.sort(
        key=lambda nitrogen: (
            mol.GetBondBetweenAtoms(guanidino_carbon.GetIdx(), nitrogen.GetIdx()).GetBondType()
            != Chem.BondType.DOUBLE,
            nitrogen.GetIdx(),
        )
    )
    return [
        (alpha_nitrogens[0].GetIdx(), "alpha_amino_N"),
        (terminal_nitrogens[0].GetIdx(), "terminal_guanidino_N1"),
        (terminal_nitrogens[1].GetIdx(), "terminal_guanidino_N2"),
    ]


def _amide_sites(partner: Partner) -> list[tuple[int, str]]:
    if partner.mol is None:
        return []
    if partner.normalized_name == "arginine":
        return _arginine_amide_sites(partner.mol)
    sites = _amine_sites(partner.mol)
    if len(sites) == 1:
        return [(sites[0], "")]
    return [(idx, f"amine_N{number}") for number, idx in enumerate(sites, start=1)]


def _alcohol_sites(mol: Chem.Mol) -> list[int]:
    sites = set()
    for match in mol.GetSubstructMatches(ALCOHOL):
        oxygen_idx = match[0]
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        carbon_neighbors = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() == 6]
        if any(_is_carboxyl_carbon(carbon) for carbon in carbon_neighbors):
            continue
        sites.add(oxygen_idx)
    return sorted(sites)


def _is_carboxyl_carbon(atom: Chem.Atom) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    has_double_oxygen = any(
        bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        for bond in atom.GetBonds()
    )
    has_single_oxygen = any(
        bond.GetBondType() == Chem.BondType.SINGLE
        and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        for bond in atom.GetBonds()
    )
    return has_double_oxygen and has_single_oxygen


def _core_hydroxyl_sites(core: BileAcidCore, requested_position: str | None = None) -> list[tuple[BileAcidSite, int]]:
    resolution = resolve_steroid_sites(core, requested_position)
    mapped_core = core.mapped_mol
    return [(site, _mapped_atom_idx(mapped_core, site.oxygen_map_num)) for site in resolution.sites]


def c24_amides(core: BileAcidCore, partner: Partner) -> list[ReactionProduct]:
    if "C24_amide" not in partner.allowed_reaction_channels or partner.mol is None:
        return []
    products = []
    mapped_core = core.mapped_mol
    carbon = _mapped_atom_idx(mapped_core, core.c24_site.carbon_map_num)
    leaving_oh = _mapped_atom_idx(mapped_core, core.c24_site.oxygen_map_num)
    acid, mapping = _remove_atom(mapped_core, leaving_oh)
    carbon = mapping[carbon]
    for n_idx, attachment_site in _amide_sites(partner):
        mol = _combine_and_bond(acid, carbon, partner.mol, n_idx)
        if mol:
            note = f"Partner attachment site: {attachment_site}." if attachment_site else ""
            if attachment_site.startswith("terminal_guanidino"):
                note += " Arginine terminal guanidino nitrogens are separate resonance/tautomer attachment models."
            products.append(
                ReactionProduct(
                    mol,
                    "C24_amide",
                    "C24_amide_BA",
                    core.c24_site,
                    "amide",
                    generation_note=note,
                    partner_attachment_site=attachment_site,
                    partner_attachment_atom_index=n_idx,
                )
            )
    return products


def c24_esters(core: BileAcidCore, partner: Partner) -> list[ReactionProduct]:
    if "C24_ester" not in partner.allowed_reaction_channels or partner.mol is None:
        return []
    products = []
    mapped_core = core.mapped_mol
    carbon = _mapped_atom_idx(mapped_core, core.c24_site.carbon_map_num)
    leaving_oh = _mapped_atom_idx(mapped_core, core.c24_site.oxygen_map_num)
    acid, mapping = _remove_atom(mapped_core, leaving_oh)
    carbon = mapping[carbon]
    for o_idx in _alcohol_sites(partner.mol)[:4]:
        mol = _combine_and_bond(acid, carbon, partner.mol, o_idx)
        if mol:
            products.append(ReactionProduct(mol, "C24_ester", "C24_ester_BA", core.c24_site, "ester"))
    return products


def ba_o_acyl(core: BileAcidCore, partner: Partner, requested_position: str | None = None) -> list[ReactionProduct]:
    if "BA_O_acyl" not in partner.allowed_reaction_channels or partner.mol is None:
        return []
    products = []
    acids = _carboxyl_sites(partner.mol)
    for site, _core_o in _core_hydroxyl_sites(core, requested_position):
        mapped_core, core_o = _core_prepared_for_oxygen_bond(core, site)
        for carbon, leaving_oh in acids[:3]:
            acyl, mapping = _remove_atom(partner.mol, leaving_oh)
            mol = _combine_and_bond(mapped_core, core_o, acyl, mapping[carbon])
            if mol:
                products.append(ReactionProduct(mol, "BA_O_acyl", "BA_O_acyl_conjugate", site, "ester"))
    return products


def _anomeric_sites(sugar: Chem.Mol) -> list[tuple[int, int]]:
    hits: list[tuple[int, int]] = []
    for atom in sugar.GetAtoms():
        if atom.GetAtomicNum() != 6 or not atom.IsInRing():
            continue
        ring_o = any(n.GetAtomicNum() == 8 and n.IsInRing() for n in atom.GetNeighbors())
        hydroxyls = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetTotalNumHs() > 0]
        if ring_o and hydroxyls:
            hits.append((atom.GetIdx(), hydroxyls[0]))
    return hits


def glycosides(core: BileAcidCore, partner: Partner, requested_position: str | None = None) -> list[ReactionProduct]:
    if "BA_O_glycoside" not in partner.allowed_reaction_channels or partner.mol is None:
        return []
    products = []
    for site, _core_o in _core_hydroxyl_sites(core, requested_position):
        mapped_core, core_o = _core_prepared_for_oxygen_bond(core, site)
        for carbon, leaving_oh in _anomeric_sites(partner.mol)[:1]:
            gly, mapping = _remove_atom(partner.mol, leaving_oh)
            mol = _combine_and_bond(mapped_core, core_o, gly, mapping[carbon])
            if mol:
                products.append(ReactionProduct(mol, "BA_O_glycoside", "BA_O_glycoside", site, "glycosidic", ["anomeric_uncertainty"]))
    return products


def glucan_fragment_products(core: BileAcidCore, max_dp: int) -> list[ReactionProduct]:
    products = []
    glucose = Chem.MolFromSmiles("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O")
    assert glucose is not None
    for dp in range(1, max_dp + 1):
        # Defined oligomer fragment model: glucose units connected by simple alpha-like O bridges.
        fragment = glucose
        for _ in range(2, dp + 1):
            fragment = _combine_and_bond(fragment, _alcohol_sites(fragment)[0], glucose, _anomeric_sites(glucose)[0][0]) or fragment
        pseudo_partner = Partner("glucan fragment", "glucan_fragment", smiles=canonical_smiles(fragment), allowed_reaction_channels=["BA_O_glycoside"])
        for product in glycosides(core, pseudo_partner):
            product.candidate_category = "BA_polymeric_fragment_model"
            product.warning_flags.extend(["polymeric_fragment_model", f"defined_degree_of_polymerization=DP{dp}", "large_flexible_warning"])
            product.generation_note = f"defined_degree_of_polymerization=DP{dp}"
            products.append(product)
    return products


def sulfate_phosphate_esters(core: BileAcidCore, requested_position: str | None = None) -> list[ReactionProduct]:
    products = []
    donors = [
        ("BA_O_sulfation", "BA_O_sulfate", "sulfate ester", "O=S(=O)(O)O", 1),
        ("BA_O_phosphorylation", "BA_O_phosphate", "phosphate ester", "O=P(O)(O)O", 1),
    ]
    for channel, category, bond_type, smiles, central in donors:
        donor = Chem.MolFromSmiles(smiles)
        if donor is None:
            continue
        leaving = [n.GetIdx() for n in donor.GetAtomWithIdx(central).GetNeighbors() if n.GetAtomicNum() == 8 and n.GetTotalNumHs() > 0][0]
        frag, mapping = _remove_atom(donor, leaving)
        for site, _core_o in _core_hydroxyl_sites(core, requested_position):
            mapped_core, core_o = _core_prepared_for_oxygen_bond(core, site)
            mol = _combine_and_bond(mapped_core, core_o, frag, mapping[central])
            if mol:
                products.append(ReactionProduct(mol, channel, category, site, bond_type))
    return products


def reactive_carbonyl_adducts(core: BileAcidCore, partner: Partner, requested_position: str | None = None) -> list[ReactionProduct]:
    if "reactive_carbonyl" not in partner.allowed_reaction_channels or partner.mol is None:
        return []
    products = []
    for c_idx, o_idx in partner.mol.GetSubstructMatches(CARBONYL)[:2]:
        rw = Chem.RWMol(partner.mol)
        rw.RemoveBond(c_idx, o_idx)
        rw.AddBond(c_idx, o_idx, Chem.BondType.SINGLE)
        hemi = rw.GetMol()
        try:
            _sanitize_and_assign(hemi)
        except Exception:
            continue
        for site, _core_o in _core_hydroxyl_sites(core, requested_position):
            mapped_core, core_o = _core_prepared_for_oxygen_bond(core, site)
            mol = _combine_and_bond(mapped_core, core_o, hemi, c_idx)
            if mol:
                products.append(
                    ReactionProduct(
                        mol,
                        "reactive_carbonyl",
                        "BA_reactive_carbonyl_adduct",
                        site,
                        "hemiacetal_or_carbonyl_adduct",
                        ["reactive_or_artifact_possible", "manual_review_required"],
                    )
                )
    return products
