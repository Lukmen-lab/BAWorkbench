from rdkit import Chem

from ba_conjugate_generator.mass import AMMONIUM_MASS, POTASSIUM_MASS, PROTON_MASS, SODIUM_MASS, calculate_mass_properties


def test_adduct_constants_are_populated():
    assert PROTON_MASS > 1.0
    assert SODIUM_MASS > 22.0
    assert AMMONIUM_MASS > 18.0
    assert POTASSIUM_MASS > 38.0


def test_formula_and_adducts_are_populated():
    props = calculate_mass_properties(Chem.MolFromSmiles("CC(=O)O"))
    assert props.formula == "C2H4O2"
    assert props.exact_mass > 60
    assert props.mz_m_plus_h > props.exact_mass
    assert props.mz_m_minus_h < props.exact_mass
