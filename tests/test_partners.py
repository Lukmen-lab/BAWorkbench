from ba_conjugate_generator.partners import built_in_partners


def by_norm():
    return {p.normalized_name: p for p in built_in_partners()}


def test_deduplicated_synonyms_are_preserved():
    partners = by_norm()
    assert "Glucose" in partners["glucose"].original_labels
    assert "D-Glucose" in partners["glucose"].original_labels
    assert "L-Glutamate" in partners["glutamic acid"].original_labels
    assert "L-Glutamic acid" in partners["glutamic acid"].original_labels
    assert "Xanthine (as artifact)" in partners["xanthine"].original_labels


def test_inorganic_salts_are_excluded_and_salt_contexts_collapsed():
    partners = by_norm()
    assert partners["acetic acid"].smiles == "CC(=O)O"
    assert "Sodium Acetate" in partners["acetic acid"].original_labels
    assert partners["pantothenic acid"].smiles
    excluded = [p.partner_name for p in built_in_partners() if p.is_excluded_context_only]
    assert any("Magnesium Sulfate" in name for name in excluded)
    assert any("Iron (II) Phosphate" in name for name in excluded)


def test_taurine_is_a_builtin_reference_c24_amide_partner():
    taurine = by_norm()["taurine"]
    assert taurine.smiles == "NCCS(=O)(=O)O"
    assert taurine.allowed_reaction_channels == ["C24_amide"]
    assert taurine.plausibility_class == "reference_known"


def test_default_partners_include_cysteine_but_not_disulfide_cystine():
    partners = by_norm()
    assert "cysteine" in partners
    assert "cystine" not in partners
