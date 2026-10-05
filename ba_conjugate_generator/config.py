"""Runtime configuration."""

from __future__ import annotations

from dataclasses import dataclass


GENERATION_MODES = {"primary_docking", "broad_covalent", "manual_review"}


@dataclass
class GeneratorConfig:
    mode: str = "broad_covalent"
    max_saccharide_dp: int = 5
    allow_dp6_to_dp10: bool = False

    def __post_init__(self) -> None:
        if self.mode not in GENERATION_MODES:
            raise ValueError(f"Unsupported generation mode: {self.mode}")
        if self.max_saccharide_dp < 1:
            raise ValueError("max_saccharide_dp must be >= 1")
        if self.max_saccharide_dp > 10 and not self.allow_dp6_to_dp10:
            raise ValueError("DP > 10 is outside the prototype model")
