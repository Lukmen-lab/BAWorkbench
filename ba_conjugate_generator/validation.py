"""Validation helpers with RDKit as the authority."""

from __future__ import annotations

from rdkit import Chem


def mol_from_smiles(smiles: str | None, name: str = "molecule") -> Chem.Mol | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def canonical_smiles(mol: Chem.Mol, isomeric: bool = True) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric, canonical=True)
