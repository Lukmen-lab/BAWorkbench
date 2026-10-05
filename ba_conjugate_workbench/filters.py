"""DataFrame filtering for the Workbench."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from .utils import split_flags


@dataclass
class FilterSpec:
    cores: list[str] = field(default_factory=list)
    categories: list[str] = field(default_factory=list)
    partners: list[str] = field(default_factory=list)
    reaction_channels: list[str] = field(default_factory=list)
    conjugation_sites: list[str] = field(default_factory=list)
    modified_positions: list[str] = field(default_factory=list)
    docking_statuses: list[str] = field(default_factory=list)
    warning_flags: list[str] = field(default_factory=list)
    mass_range: tuple[float | None, float | None] = (None, None)
    mz_minus_h_range: tuple[float | None, float | None] = (None, None)
    mz_plus_h_range: tuple[float | None, float | None] = (None, None)
    formula_contains_n: bool = False
    formula_contains_s: bool = False
    formula_contains_p: bool = False
    manual_review_required: bool | None = None
    text_search: str = ""


def _isin(df: pd.DataFrame, column: str, values: list[str]) -> pd.Series:
    if not values or column not in df.columns:
        return pd.Series(True, index=df.index)
    return df[column].astype(str).isin(values)


def _range(df: pd.DataFrame, column: str, bounds: tuple[float | None, float | None]) -> pd.Series:
    if column not in df.columns:
        return pd.Series(True, index=df.index)
    values = pd.to_numeric(df[column], errors="coerce")
    mask = pd.Series(True, index=df.index)
    low, high = bounds
    if low is not None:
        mask &= values >= low
    if high is not None:
        mask &= values <= high
    return mask


def apply_filters(df: pd.DataFrame, spec: FilterSpec) -> pd.DataFrame:
    mask = pd.Series(True, index=df.index)
    mask &= _isin(df, "parent_core_abbrev", spec.cores)
    mask &= _isin(df, "candidate_category", spec.categories)
    mask &= _isin(df, "partner_name", spec.partners)
    mask &= _isin(df, "reaction_channel", spec.reaction_channels)
    mask &= _isin(df, "conjugation_site", spec.conjugation_sites)
    mask &= _isin(df, "modified_position", spec.modified_positions)
    mask &= _isin(df, "docking_status", spec.docking_statuses)
    mask &= _range(df, "neutral_monoisotopic_mass", spec.mass_range)
    mask &= _range(df, "mz_M_minus_H", spec.mz_minus_h_range)
    mask &= _range(df, "mz_M_plus_H", spec.mz_plus_h_range)
    if spec.warning_flags and "warning_flags" in df.columns:
        mask &= df["warning_flags"].apply(lambda v: any(flag in split_flags(v) for flag in spec.warning_flags))
    if spec.formula_contains_n and "molecular_formula" in df.columns:
        mask &= df["molecular_formula"].astype(str).str.contains("N", regex=False)
    if spec.formula_contains_s and "molecular_formula" in df.columns:
        mask &= df["molecular_formula"].astype(str).str.contains("S", regex=False)
    if spec.formula_contains_p and "molecular_formula" in df.columns:
        mask &= df["molecular_formula"].astype(str).str.contains("P", regex=False)
    if spec.manual_review_required is not None and "manual_review_required" in df.columns:
        wanted = str(spec.manual_review_required).lower()
        mask &= df["manual_review_required"].astype(str).str.lower().isin({wanted, "1" if spec.manual_review_required else "0"})
    if spec.text_search:
        columns = [c for c in ["candidate_id", "candidate_name", "partner_name", "isomeric_smiles", "canonical_smiles", "MS_notes", "docking_notes"] if c in df.columns]
        haystack = df[columns].astype(str).agg(" ".join, axis=1).str.lower()
        mask &= haystack.str.contains(spec.text_search.lower(), regex=False)
    return df.loc[mask].copy()


def column_options(df: pd.DataFrame, column: str) -> list[str]:
    if column not in df.columns:
        return []
    return sorted(v for v in df[column].astype(str).unique().tolist() if v and v.lower() != "nan")


def warning_flag_options(df: pd.DataFrame) -> list[str]:
    if "warning_flags" not in df.columns:
        return []
    flags: set[str] = set()
    for value in df["warning_flags"]:
        flags.update(split_flags(value))
    return sorted(flags)
