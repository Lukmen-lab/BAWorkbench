"""QC summary helpers."""

from __future__ import annotations

from collections import Counter

import pandas as pd

from .utils import split_flags


def qc_summary(qc: dict[str, object], candidates: pd.DataFrame | None = None) -> dict[str, object]:
    summary: dict[str, object] = {}
    for key in ["invalid_candidates", "excluded_partners", "manual_review_candidates", "duplicate_or_isomeric_groups"]:
        value = qc.get(key)
        summary[f"{key}_count"] = len(value) if isinstance(value, pd.DataFrame) else None
    if candidates is not None and not candidates.empty:
        flags = Counter()
        if "warning_flags" in candidates.columns:
            for value in candidates["warning_flags"]:
                flags.update(split_flags(value))
        summary["most_common_warning_flags"] = flags.most_common(10)
        if "candidate_category" in candidates.columns:
            summary["most_common_candidate_categories"] = candidates["candidate_category"].value_counts().head(10).to_dict()
    gen_summary = qc.get("generation_summary")
    if isinstance(gen_summary, dict):
        summary["generation_summary"] = gen_summary
    return summary
