"""Candidate and QC file loading."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .utils import REQUIRED_COLUMNS, SMILES_COLUMNS, first_existing


CANDIDATE_FILENAMES = (
    "candidates_full.tsv",
    "candidates_full.csv",
    "candidates_metadata.json",
)


@dataclass
class LoadResult:
    candidates: pd.DataFrame
    source_path: Path
    warnings: list[str] = field(default_factory=list)


def _normalize_candidate_frame(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    missing_required = REQUIRED_COLUMNS - set(df.columns)
    if missing_required:
        raise ValueError(f"Candidate table is missing required columns: {', '.join(sorted(missing_required))}")
    if not any(col in df.columns for col in SMILES_COLUMNS):
        raise ValueError("Candidate table requires either isomeric_smiles or canonical_smiles")
    for column in SMILES_COLUMNS:
        if column not in df.columns:
            warnings.append(f"Optional SMILES column missing: {column}")
            df[column] = ""
    if "curation_status" not in df.columns:
        df["curation_status"] = "unreviewed"
    if "curation_notes" not in df.columns:
        df["curation_notes"] = ""
    for column in ["neutral_monoisotopic_mass", "mz_M_minus_H", "mz_M_plus_H"]:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df.fillna(""), warnings


def load_candidate_file(path: str | Path) -> LoadResult:
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    suffix = source.suffix.lower()
    if suffix == ".json":
        records = json.loads(source.read_text(encoding="utf-8"))
        df = pd.DataFrame.from_records(records)
    elif suffix == ".csv":
        df = pd.read_csv(source)
    else:
        df = pd.read_csv(source, sep="\t")
    df, warnings = _normalize_candidate_frame(df)
    return LoadResult(df, source, warnings)


def load_run_folder(path: str | Path) -> LoadResult:
    folder = Path(path)
    source = first_existing([folder / name for name in CANDIDATE_FILENAMES])
    if source is None:
        raise FileNotFoundError(f"No candidates_full.tsv/csv/json found in {folder}")
    return load_candidate_file(source)


def discover_run_folders(search_root: str | Path = "outputs") -> list[str]:
    """Return immediate child folders that contain a loadable candidate table."""
    root = Path(search_root)
    if not root.is_dir():
        return []
    folders = [
        folder.as_posix()
        for folder in root.iterdir()
        if folder.is_dir() and any((folder / name).is_file() for name in CANDIDATE_FILENAMES)
    ]
    return sorted(folders, key=str.casefold)


def load_qc_files(run_folder: str | Path) -> dict[str, object]:
    folder = Path(run_folder)
    qc: dict[str, object] = {}
    table_files = {
        "invalid_candidates": "invalid_candidates.tsv",
        "excluded_partners": "excluded_partners.tsv",
        "manual_review_candidates": "manual_review_candidates.tsv",
        "duplicate_or_isomeric_groups": "duplicate_or_isomeric_groups.tsv",
    }
    for key, name in table_files.items():
        path = folder / name
        qc[key] = pd.read_csv(path, sep="\t") if path.exists() else None
    summary_path = folder / "generation_summary.json"
    qc["generation_summary"] = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else None
    return qc
