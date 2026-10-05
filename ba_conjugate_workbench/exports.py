"""Workbench export helpers."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from rdkit import Chem

from .conformer_cache import generate_conformer
from .utils import candidate_smiles, ensure_directory


def export_table(df: pd.DataFrame, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if out.suffix.lower() == ".csv" else "\t"
    df.to_csv(out, sep=sep, index=False)
    return out


def export_metadata_json(df: pd.DataFrame, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(df.to_dict(orient="records"), indent=2), encoding="utf-8")
    return out


def export_docking_smiles(df: pd.DataFrame, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, row in df.iterrows():
        rows.append(
            {
                "candidate_id": row.get("candidate_id", ""),
                "isomeric_smiles": row.get("isomeric_smiles") or row.get("canonical_smiles", ""),
                "candidate_name": row.get("candidate_name", row.get("candidate_id", "")),
            }
        )
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return out


def export_manual_review_notes(df: pd.DataFrame, path: str | Path) -> Path:
    if "manual_review_required" in df.columns:
        mask = df["manual_review_required"].astype(str).str.lower().isin({"true", "1", "yes"})
        out_df = df.loc[mask].copy()
    else:
        out_df = df[df.get("curation_status", "") == "manual_review"].copy()
    return export_table(out_df, path)


def export_sdf(df: pd.DataFrame, path: str | Path, cache_dir: str | Path) -> tuple[Path, list[str]]:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    messages: list[str] = []
    writer = Chem.SDWriter(str(out))
    try:
        for _, row in df.iterrows():
            smiles = candidate_smiles(row.to_dict())
            if not smiles:
                messages.append(f"Skipped {row.get('candidate_id', '')}: missing SMILES")
                continue
            result = generate_conformer(smiles, cache_dir=cache_dir)
            if result.mol is None:
                messages.append(f"Skipped {row.get('candidate_id', '')}: {result.error}")
                continue
            result.mol.SetProp("_Name", str(row.get("candidate_id", "")))
            for col in ["candidate_name", "molecular_formula", "neutral_monoisotopic_mass", "docking_status"]:
                if col in row:
                    result.mol.SetProp(col, str(row[col]))
            writer.write(result.mol)
    finally:
        writer.close()
    return out, messages


def export_workbench_bundle(visible: pd.DataFrame, selected: pd.DataFrame, curated: pd.DataFrame, out_dir: str | Path) -> dict[str, Path]:
    folder = ensure_directory(out_dir)
    return {
        "visible_candidates": export_table(visible, folder / "visible_candidates.tsv"),
        "selected_candidates": export_table(selected, folder / "selected_candidates.tsv"),
        "curated_candidates_full": export_table(curated, folder / "curated_candidates_full.tsv"),
        "curated_smiles_for_docking": export_docking_smiles(curated, folder / "curated_smiles_for_docking.tsv"),
        "curated_candidates_metadata": export_metadata_json(curated, folder / "curated_candidates_metadata.json"),
        "manual_review_notes": export_manual_review_notes(curated, folder / "manual_review_notes.tsv"),
    }
