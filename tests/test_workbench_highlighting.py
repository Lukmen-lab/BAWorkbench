from ba_conjugate_workbench.data_loader import load_candidate_file
from ba_conjugate_workbench.highlighting import conventional_ba_labels, infer_highlights, large_flexible_warning


def test_highlighting_defaults_when_metadata_missing():
    spec = infer_highlights({"candidate_id": "x", "isomeric_smiles": "CCO"})
    assert isinstance(spec.core_atoms, list)
    assert isinstance(spec.warnings, list)


def test_highlighting_for_generated_candidate():
    df = load_candidate_file("outputs/run_001/candidates_full.tsv").candidates
    row = df.iloc[0].to_dict()
    spec = infer_highlights(row)
    assert spec.core_atoms or spec.site_atoms


def test_large_flexible_warning_flag():
    warning = large_flexible_warning({"warning_flags": "large_flexible_warning", "neutral_monoisotopic_mass": 1000})
    assert warning is not None


def test_conventional_ba_labels_are_not_rdkit_position_numbers():
    df = load_candidate_file("outputs/run_001/candidates_full.tsv").candidates
    row = df[df["candidate_id"] == "CA_0b261eba308e"].iloc[0].to_dict()
    labels = conventional_ba_labels(row)
    assert "C24" in labels.values()
    assert {"C3", "C7", "C12", "C24"}.issubset(set(labels.values()))
    assert all(not label.startswith("C0") for label in labels.values())
