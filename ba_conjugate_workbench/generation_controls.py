"""Thin GUI-facing wrapper around the existing generator."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ba_conjugate_generator.config import GeneratorConfig
from ba_conjugate_generator.cores import get_cores
from ba_conjugate_generator.exports import export_outputs
from ba_conjugate_generator.generator import generate_candidates, unmodified_record
from ba_conjugate_generator.grouping import assign_structure_duplicate_groups, group_ids
from ba_conjugate_generator.partners import load_partners


@dataclass
class WorkbenchGenerationRequest:
    cores: list[str]
    mode: str
    output: str | Path
    max_saccharide_dp: int = 5
    partners_path: str | Path | None = None
    include_reactive: bool = True
    include_large_flexible: bool = True
    include_manual_review: bool = True


def run_generation(request: WorkbenchGenerationRequest) -> dict[str, object]:
    config = GeneratorConfig(mode=request.mode, max_saccharide_dp=request.max_saccharide_dp, allow_dp6_to_dp10=True)
    cores = get_cores(request.cores)
    partners = load_partners(request.partners_path)
    if not request.include_reactive:
        partners = [p for p in partners if "reactive_carbonyl" not in p.allowed_reaction_channels]
    candidates, invalid, excluded = generate_candidates(cores, partners, config)
    if not request.include_large_flexible:
        candidates = [c for c in candidates if "large_flexible_warning" not in str(c.get("warning_flags", ""))]
    if not request.include_manual_review:
        candidates = [c for c in candidates if not c.get("manual_review_required")]
    repaired_unmodified = False
    for core in cores:
        matching = [
            candidate
            for candidate in candidates
            if candidate.get("reaction_channel") == "unmodified"
            and candidate.get("parent_core_abbrev") == core.abbrev
        ]
        if not matching:
            candidates.append(unmodified_record(core, config))
            repaired_unmodified = True
    if repaired_unmodified:
        group_ids(candidates)
        assign_structure_duplicate_groups(candidates)
    unmodified = [c for c in candidates if c.get("reaction_channel") == "unmodified"]
    unmodified_counts = {
        core.abbrev: sum(c.get("parent_core_abbrev") == core.abbrev for c in unmodified)
        for core in cores
    }
    export_outputs(request.output, candidates, invalid, excluded)
    return {
        "output": str(request.output),
        "candidate_count": len(candidates),
        "unmodified_ba_count": len(unmodified),
        "unmodified_ba_repaired": repaired_unmodified,
        "invalid_count": len(invalid),
        "excluded_count": len(excluded),
    }
