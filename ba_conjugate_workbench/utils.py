"""Shared Workbench utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIR = PROJECT_ROOT / "outputs" / "run_001"
DEFAULT_CACHE_DIR = PROJECT_ROOT / "cache" / "conformer_cache"

REQUIRED_COLUMNS = {"candidate_id"}
SMILES_COLUMNS = ("isomeric_smiles", "canonical_smiles")


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def candidate_smiles(row: dict) -> str:
    for column in SMILES_COLUMNS:
        value = row.get(column)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def split_flags(value: object) -> list[str]:
    if value is None:
        return []
    text = str(value)
    if not text or text.lower() == "nan":
        return []
    return [part.strip() for part in text.replace(",", ";").split(";") if part.strip()]


def ensure_directory(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out
