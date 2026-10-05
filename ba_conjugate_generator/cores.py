"""Bile acid core registry."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from .validation import mol_from_smiles


@dataclass(frozen=True)
class BileAcidSite:
    position: str
    conjugation_site: str
    carbon_map_num: int
    oxygen_map_num: int
    site_kind: str


SITE_DEFINITIONS = {
    "C3": BileAcidSite("C3", "C3-OH", 3, 103, "steroid_hydroxyl"),
    "C7": BileAcidSite("C7", "C7-OH", 7, 107, "steroid_hydroxyl"),
    "C12": BileAcidSite("C12", "C12-OH", 12, 112, "steroid_hydroxyl"),
    "C24": BileAcidSite("C24", "C24-COOH", 24, 124, "sidechain_carboxyl"),
}


@dataclass(frozen=True)
class BileAcidCore:
    abbrev: str
    name: str
    smiles: str
    sites: tuple[str, ...]
    hydroxyl_sites: tuple[str, ...]
    notes: str = ""
    mapped_smiles: str = ""

    @property
    def mol(self) -> Chem.Mol:
        mol = mol_from_smiles(self.smiles, self.abbrev)
        if mol is None:
            raise ValueError(f"Invalid core SMILES for {self.abbrev}")
        return mol

    @property
    def mapped_mol(self) -> Chem.Mol:
        mol = mol_from_smiles(self.mapped_smiles, f"{self.abbrev} mapped template")
        if mol is None:
            raise ValueError(f"Invalid mapped core SMILES for {self.abbrev}")
        return mol

    @property
    def hydroxyl_site_objects(self) -> tuple[BileAcidSite, ...]:
        positions = tuple(name.split("-")[0] for name in self.hydroxyl_sites)
        return tuple(SITE_DEFINITIONS[position] for position in positions)

    @property
    def c24_site(self) -> BileAcidSite:
        return SITE_DEFINITIONS["C24"]

    def site(self, position: str) -> BileAcidSite | None:
        normalized = position.replace("-OH", "").replace("-COOH", "")
        if normalized == "C24":
            return self.c24_site
        for site in self.hydroxyl_site_objects:
            if site.position == normalized:
                return site
        return None

    def atom_index(self, atom_map_num: int) -> int:
        matches = [atom.GetIdx() for atom in self.mapped_mol.GetAtoms() if atom.GetAtomMapNum() == atom_map_num]
        if len(matches) != 1:
            raise ValueError(f"{self.abbrev} template has {len(matches)} atoms with map {atom_map_num}")
        return matches[0]


CORE_REGISTRY: dict[str, BileAcidCore] = {
    "CA": BileAcidCore(
        "CA",
        "cholic acid",
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
        ("C3-OH", "C7-OH", "C12-OH", "C24-COOH"),
        ("C3-OH", "C7-OH", "C12-OH"),
        "SMILES seeded from public structure records; validate before production use.",
        "[CH2:1]1[CH2:2][C@@H:3]([OH:103])[CH2:4][C@H:5]2[CH2:6][C@@H:7]([OH:107])[C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][C@H:12]([OH:112])[C@@:13]1([CH3:18])[C@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
    "CDCA": BileAcidCore(
        "CDCA",
        "chenodeoxycholic acid",
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C",
        ("C3-OH", "C7-OH", "C24-COOH"),
        ("C3-OH", "C7-OH"),
        "SMILES seeded from public structure records; validate before production use.",
        "[CH2:1]1[CH2:2][C@@H:3]([OH:103])[CH2:4][C@H:5]2[CH2:6][C@@H:7]([OH:107])[C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][CH2:12][C@@:13]1([CH3:18])[C@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
    "DCA": BileAcidCore(
        "DCA",
        "deoxycholic acid",
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C",
        ("C3-OH", "C12-OH", "C24-COOH"),
        ("C3-OH", "C12-OH"),
        "SMILES seeded from public structure records; validate before production use.",
        "[CH2:1]1[CH2:2][C@@H:3]([OH:103])[CH2:4][C@H:5]2[CH2:6][CH2:7][C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][C@H:12]([OH:112])[C@@:13]1([CH3:18])[C@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
    "UDCA": BileAcidCore(
        "UDCA",
        "ursodeoxycholic acid",
        "O=C(O)CC[C@H]([C@H]1CC[C@@H]2[C@]1(C)CC[C@H]4[C@H]2[C@@H](O)C[C@@H]3C[C@H](O)CC[C@@]34C)C",
        ("C3-OH", "C7-OH", "C24-COOH"),
        ("C3-OH", "C7-OH"),
        "Preserves 7-hydroxyl stereochemistry distinct from CDCA.",
        "[CH2:1]1[CH2:2][C@@H:3]([OH:103])[CH2:4][C@H:5]2[CH2:6][C@H:7]([OH:107])[C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][CH2:12][C@@:13]1([CH3:18])[C@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
    "LCA": BileAcidCore(
        "LCA",
        "lithocholic acid",
        "O=C(O)CC[C@H]([C@H]1CC[C@@H]2[C@]1(C)CC[C@H]4[C@H]2CC[C@@H]3C[C@H](O)CC[C@@]34C)C",
        ("C3-OH", "C24-COOH"),
        ("C3-OH",),
        "3alpha-hydroxy-5beta-cholan-24-oic acid; stereochemistry checked against PubChem CID 9903.",
        "[CH2:1]1[CH2:2][C@@H:3]([OH:103])[CH2:4][C@H:5]2[CH2:6][CH2:7][C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][CH2:12][C@@:13]1([CH3:18])[C@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
    "7oxoCA": BileAcidCore(
        "7oxoCA",
        "7-oxo-cholic acid",
        "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2C(=O)C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C",
        ("C3-OH", "C7-oxo", "C12-OH", "C24-COOH"),
        ("C3-OH", "C12-OH"),
        "TODO: confirm exact stereochemical source structure before production use.",
        "[CH2:1]1[CH2:2][C@@H:3]([OH:103])[CH2:4][C@H:5]2[CH2:6][C:7](=[O:207])[C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][C@H:12]([OH:112])[C@@:13]1([CH3:18])[C@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
    "7oxoCDCA": BileAcidCore(
        "7oxoCDCA",
        "7-oxo-chenodeoxycholic acid",
        "O=C(O)CCC(C)C1CCC2C3C(=O)CC4CC(O)CCC4(C)C3CCC12C",
        ("C3-OH", "C7-oxo", "C24-COOH"),
        ("C3-OH",),
        "TODO: stereochemistry incomplete in source SMILES; flagged for review.",
        "[CH2:1]1[CH2:2][CH:3]([OH:103])[CH2:4][CH:5]2[CH2:6][C:7](=[O:207])[CH:8]3[CH:9]([C:10]12[CH3:19])[CH2:11][CH2:12][C:13]1([CH3:18])[CH:14]3[CH2:15][CH2:16][CH:17]1[CH:20]([CH3:21])[CH2:22][CH2:23][C:24]([OH:124])=[O:224]",
    ),
}


def get_cores(selected: list[str] | None = None) -> list[BileAcidCore]:
    keys = selected or list(CORE_REGISTRY)
    missing = [key for key in keys if key not in CORE_REGISTRY]
    if missing:
        raise KeyError(f"Unknown core(s): {', '.join(missing)}")
    return [CORE_REGISTRY[key] for key in keys]
