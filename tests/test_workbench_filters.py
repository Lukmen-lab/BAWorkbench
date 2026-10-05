from ba_conjugate_workbench.data_loader import load_candidate_file
from ba_conjugate_workbench.filters import FilterSpec, apply_filters


def test_filter_by_core_category_reaction_channel():
    df = load_candidate_file("outputs/run_001/candidates_full.tsv").candidates
    spec = FilterSpec(cores=["CA"], reaction_channels=["C24_amide"])
    filtered = apply_filters(df, spec)
    assert not filtered.empty
    assert set(filtered["parent_core_abbrev"]) == {"CA"}
    assert set(filtered["reaction_channel"]) == {"C24_amide"}


def test_text_and_formula_filters():
    df = load_candidate_file("outputs/run_001/candidates_full.tsv").candidates
    filtered = apply_filters(df, FilterSpec(text_search="glycine", formula_contains_n=True))
    assert not filtered.empty
    assert filtered["molecular_formula"].str.contains("N").all()
