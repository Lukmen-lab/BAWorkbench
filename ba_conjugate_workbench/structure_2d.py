"""RDKit 2D depiction helpers."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

from .highlighting import HighlightSpec, infer_highlights
from .utils import candidate_smiles


def mol_from_candidate(row: dict) -> Chem.Mol | None:
    smiles = candidate_smiles(row)
    return Chem.MolFromSmiles(smiles) if smiles else None


def atom_label(atom: Chem.Atom, mode: str, conventional: str | None = None) -> str:
    if mode == "Element only":
        return atom.GetSymbol()
    if mode == "Element + RDKit atom index":
        return f"{atom.GetSymbol()}{atom.GetIdx()}"
    if mode == "BA conventional labels + atom indices":
        return f"{conventional} (idx {atom.GetIdx()})" if conventional else ""
    return ""


def depiction_svg(row: dict, label_mode: str = "Off", width: int = 700, height: int = 420) -> tuple[str | None, list[str]]:
    mol = mol_from_candidate(row)
    if mol is None:
        return None, ["Invalid or missing SMILES; 2D depiction unavailable."]
    mol = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol)
    spec = infer_highlights(row)
    if label_mode != "Off":
        for atom in mol.GetAtoms():
            conventional = spec.labels.get(atom.GetIdx())
            label = atom_label(atom, label_mode, conventional)
            if label:
                atom.SetProp("atomLabel", label)
    drawer = Draw.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.padding = 0.08
    highlight_atoms = sorted(set(spec.core_atoms + spec.partner_atoms + spec.site_atoms))
    atom_colours = {idx: (0.35, 0.55, 0.95) for idx in spec.core_atoms}
    atom_colours.update({idx: (0.95, 0.60, 0.20) for idx in spec.partner_atoms})
    atom_colours.update({idx: (0.95, 0.20, 0.25) for idx in spec.site_atoms})
    bond_colours = {idx: (0.95, 0.10, 0.10) for idx in spec.new_bonds}
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=spec.new_bonds,
        highlightAtomColors=atom_colours,
        highlightBondColors=bond_colours,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText(), spec.warnings
