"""On-demand RDKit conformer generation and disk cache."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import DEFAULT_CACHE_DIR


@dataclass
class ConformerResult:
    mol: Chem.Mol | None
    mol_block: str
    cache_path: Path | None
    from_cache: bool
    status_messages: list[str]
    error: str | None = None


def cache_key(smiles: str, show_hydrogens: bool = True, method: str = "MMFF", seed: int = 61453) -> str:
    raw = f"{smiles}|H={show_hydrogens}|method={method}|seed={seed}"
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]


def cache_path_for(smiles: str, cache_dir: str | Path = DEFAULT_CACHE_DIR, show_hydrogens: bool = True, method: str = "MMFF", seed: int = 61453) -> Path:
    return Path(cache_dir) / f"{cache_key(smiles, show_hydrogens, method, seed)}.sdf"


def load_cached_conformer(path: str | Path) -> Chem.Mol | None:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    return supplier[0] if supplier and len(supplier) and supplier[0] is not None else None


def save_conformer(mol: Chem.Mol, path: str | Path) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out))
    writer.write(mol)
    writer.close()


def generate_conformer(
    smiles: str,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    show_hydrogens: bool = True,
    method: str = "MMFF",
    seed: int = 61453,
    use_cache: bool = True,
) -> ConformerResult:
    path = cache_path_for(smiles, cache_dir, show_hydrogens, method, seed)
    messages: list[str] = []
    if use_cache and path.exists():
        mol = load_cached_conformer(path)
        if mol is not None:
            messages.append("Cached conformer loaded.")
            return ConformerResult(mol, Chem.MolToMolBlock(mol), path, True, messages)
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        return ConformerResult(None, "", None, False, messages, "3D generation failed: invalid SMILES.")
    mol = Chem.AddHs(base)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return ConformerResult(None, "", None, False, messages, "3D generation failed: RDKit embedding did not converge.")
    messages.append("3D conformer generated successfully.")
    optimized = False
    if method.upper() == "MMFF" and AllChem.MMFFHasAllMoleculeParams(mol):
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
            messages.append("MMFF optimization successful.")
            optimized = True
        except Exception:
            messages.append("MMFF optimization failed; UFF fallback attempted.")
    else:
        messages.append("MMFF unavailable; UFF used.")
    if not optimized:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=300)
            messages.append("UFF optimization successful.")
        except Exception as exc:
            messages.append(f"UFF optimization failed: {exc}")
    output_mol = mol if show_hydrogens else Chem.RemoveHs(mol)
    save_conformer(output_mol, path)
    return ConformerResult(output_mol, Chem.MolToMolBlock(output_mol), path, False, messages)


def clear_cache(cache_dir: str | Path = DEFAULT_CACHE_DIR) -> None:
    path = Path(cache_dir)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def cache_size(cache_dir: str | Path = DEFAULT_CACHE_DIR) -> tuple[int, int]:
    path = Path(cache_dir)
    if not path.exists():
        return 0, 0
    files = list(path.glob("*.sdf"))
    return len(files), sum(file.stat().st_size for file in files)
