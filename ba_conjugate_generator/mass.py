"""Mass and adduct calculations."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


PROTON_MASS = 1.007276466621
SODIUM_MASS = 22.989218
AMMONIUM_MASS = 18.033823
POTASSIUM_MASS = 38.963158


@dataclass(frozen=True)
class MassProperties:
    formula: str
    exact_mass: float
    formal_charge: int
    mz_m_minus_h: float
    mz_m_plus_h: float
    mz_m_plus_na: float
    mz_m_plus_nh4: float
    mz_m_plus_k: float


def calculate_mass_properties(mol: Chem.Mol) -> MassProperties:
    formula = rdMolDescriptors.CalcMolFormula(mol)
    exact_mass = Descriptors.ExactMolWt(mol)
    charge = Chem.GetFormalCharge(mol)
    return MassProperties(
        formula=formula,
        exact_mass=exact_mass,
        formal_charge=charge,
        mz_m_minus_h=exact_mass - PROTON_MASS,
        mz_m_plus_h=exact_mass + PROTON_MASS,
        mz_m_plus_na=exact_mass + SODIUM_MASS,
        mz_m_plus_nh4=exact_mass + AMMONIUM_MASS,
        mz_m_plus_k=exact_mass + POTASSIUM_MASS,
    )
