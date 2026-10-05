from ba_conjugate_generator.config import GeneratorConfig
from ba_conjugate_generator.cores import get_cores
from ba_conjugate_generator.generator import generate_candidates
from ba_conjugate_generator.partners import built_in_partners
from ba_conjugate_generator.reactions import _alcohol_sites


def _partner(name):
    return next(p for p in built_in_partners() if p.normalized_name == name)


def test_unmodified_parent_bile_acids_are_included_without_partners():
    rows, invalid, excluded = generate_candidates(
        get_cores(["CA", "LCA"]),
        [],
        GeneratorConfig(mode="primary_docking"),
    )

    assert not invalid
    assert not excluded
    assert len(rows) == 2
    assert {row["parent_core_abbrev"] for row in rows} == {"CA", "LCA"}
    assert {row["molecular_formula"] for row in rows} == {"C24H40O5", "C24H40O3"}
    assert all(row["candidate_category"] == "unmodified_BA" for row in rows)
    assert all(row["reaction_channel"] == "unmodified" for row in rows)
    assert all(row["modified_position"] == "none" for row in rows)
    assert all(row["mass_shift_from_parent_core"] == 0.0 for row in rows)
    assert all(row["docking_ready"] for row in rows)


def test_c24_glycine_amide_generation_works_for_ca():
    rows, invalid, _ = generate_candidates(get_cores(["CA"]), [_partner("glycine")], GeneratorConfig(mode="primary_docking"))
    assert not invalid
    assert any(r["reaction_channel"] == "C24_amide" and r["molecular_formula"] for r in rows)


def test_builtin_taurine_generates_taurocholic_acid():
    rows, invalid, _ = generate_candidates(get_cores(["CA"]), [_partner("taurine")], GeneratorConfig(mode="primary_docking"))
    assert not invalid
    tauroca = next(r for r in rows if r["partner_normalized_name"] == "taurine" and r["reaction_channel"] == "C24_amide")
    assert tauroca["molecular_formula"] == "C26H45NO7S"
    assert tauroca["modified_position"] == "C24"
    assert tauroca["site_validation_passed"]


def test_builtin_taurine_generates_taurolithocholic_acid():
    rows, invalid, _ = generate_candidates(get_cores(["LCA"]), [_partner("taurine")], GeneratorConfig(mode="primary_docking"))
    assert not invalid
    taurolca = next(r for r in rows if r["partner_normalized_name"] == "taurine" and r["reaction_channel"] == "C24_amide")
    assert taurolca["molecular_formula"] == "C26H45NO5S"
    assert taurolca["modified_position"] == "C24"
    assert taurolca["site_validation_passed"]


def test_c24_l_phenylalanine_amide_generation_works():
    rows, _, _ = generate_candidates(get_cores(["CA"]), [_partner("phenylalanine")], GeneratorConfig(mode="primary_docking"))
    assert any(r["reaction_channel"] == "C24_amide" for r in rows)


def test_arginine_amides_use_alpha_and_both_terminal_guanidino_nitrogens():
    rows, invalid, _ = generate_candidates(
        get_cores(["CA"]),
        [_partner("arginine")],
        GeneratorConfig(mode="primary_docking"),
    )
    arginine_amides = [
        row
        for row in rows
        if row["partner_normalized_name"] == "arginine" and row["reaction_channel"] == "C24_amide"
    ]

    assert not invalid
    assert len(arginine_amides) == 3
    assert {row["partner_attachment_site"] for row in arginine_amides} == {
        "alpha_amino_N",
        "terminal_guanidino_N1",
        "terminal_guanidino_N2",
    }
    assert len({row["candidate_name"] for row in arginine_amides}) == 3
    assert len({row["isomeric_smiles"] for row in arginine_amides}) == 3
    assert all(row["partner_attachment_site"] in row["candidate_name"] for row in arginine_amides)


def test_ba_o_acetylation_works_using_acetate():
    rows, _, _ = generate_candidates(get_cores(["CA"]), [_partner("acetic acid")], GeneratorConfig(mode="primary_docking"))
    assert any(r["reaction_channel"] == "BA_O_acyl" and r["bond_type"] == "ester" for r in rows)


def test_ba_o_glycoside_generation_works_with_d_glucose():
    rows, _, _ = generate_candidates(get_cores(["CA"]), [_partner("glucose")], GeneratorConfig(mode="primary_docking"))
    assert any(r["reaction_channel"] == "BA_O_glycoside" for r in rows)


def test_disaccharides_and_glucan_fragments_in_broad_mode():
    ps = [_partner("lactose"), _partner("maltose"), _partner("sucrose"), _partner("glucan_exopolysaccharide")]
    rows, _, _ = generate_candidates(get_cores(["CA"]), ps, GeneratorConfig(mode="broad_covalent", max_saccharide_dp=3))
    assert any(r["partner_normalized_name"] == "lactose" for r in rows)
    assert any(r["partner_normalized_name"] == "maltose" for r in rows)
    assert any(r["partner_normalized_name"] == "sucrose" for r in rows)
    assert any("defined_degree_of_polymerization=DP3" in r["warning_flags"] for r in rows)


def test_reactive_candidates_are_flagged():
    rows, _, _ = generate_candidates(get_cores(["CA"]), [_partner("acrolein")], GeneratorConfig(mode="broad_covalent"))
    reactive = [r for r in rows if r["partner_normalized_name"] == "acrolein"]
    assert reactive
    assert all(r["manual_review_required"] for r in reactive)
    assert all(r["docking_status"] == "manual_review_required" for r in reactive)


def test_core_alcohol_sites_exclude_c24_carboxylic_acid_oh():
    ca = get_cores(["CA"])[0]
    assert len(_alcohol_sites(ca.mol)) == 3
    assert len(_alcohol_sites(get_cores(["LCA"])[0].mol)) == 1
