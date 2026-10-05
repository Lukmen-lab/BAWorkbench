"""Curation state management."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


CURATION_STATUSES = [
    "unreviewed",
    "include_for_docking",
    "exclude_from_docking",
    "manual_review",
    "high_priority_for_MS",
    "interesting_but_later",
]


CURATION_COLUMNS = ["candidate_id", "curation_status", "curation_notes", "reviewer", "review_timestamp"]


def empty_curation(candidate_ids: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": candidate_ids,
            "curation_status": "unreviewed",
            "curation_notes": "",
            "reviewer": "",
            "review_timestamp": "",
        }
    )


def load_curation(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=CURATION_COLUMNS)
    df = pd.read_csv(p, sep="\t")
    for col in CURATION_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[CURATION_COLUMNS].fillna("")


def save_curation(df: pd.DataFrame, path: str | Path) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df[CURATION_COLUMNS].to_csv(out, sep="\t", index=False)


def update_curation(df: pd.DataFrame, candidate_id: str, status: str, notes: str, reviewer: str) -> pd.DataFrame:
    if status not in CURATION_STATUSES:
        raise ValueError(f"Unsupported curation status: {status}")
    if df.empty or "candidate_id" not in df.columns:
        df = pd.DataFrame(columns=CURATION_COLUMNS)
    df = df.copy()
    timestamp = datetime.now(timezone.utc).isoformat()
    mask = df["candidate_id"].astype(str) == str(candidate_id)
    values = {
        "candidate_id": candidate_id,
        "curation_status": status,
        "curation_notes": notes,
        "reviewer": reviewer,
        "review_timestamp": timestamp,
    }
    if mask.any():
        for key, value in values.items():
            df.loc[mask, key] = value
    else:
        df = pd.concat([df, pd.DataFrame([values])], ignore_index=True)
    return df[CURATION_COLUMNS]


def merge_curation(candidates: pd.DataFrame, curation: pd.DataFrame) -> pd.DataFrame:
    base = candidates.drop(columns=[c for c in ["curation_status", "curation_notes", "reviewer", "review_timestamp"] if c in candidates.columns])
    if curation.empty:
        curation = empty_curation(base["candidate_id"].astype(str).tolist())
    merged = base.merge(curation, how="left", on="candidate_id")
    merged["curation_status"] = merged["curation_status"].fillna("unreviewed")
    merged["curation_notes"] = merged["curation_notes"].fillna("")
    merged["reviewer"] = merged["reviewer"].fillna("")
    merged["review_timestamp"] = merged["review_timestamp"].fillna("")
    return merged
