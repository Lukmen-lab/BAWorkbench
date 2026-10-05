import pandas as pd

from ba_conjugate_workbench.data_loader import discover_run_folders, load_candidate_file, load_qc_files, load_run_folder


def test_load_candidates_full_tsv():
    result = load_candidate_file("outputs/run_001/candidates_full.tsv")
    assert len(result.candidates) > 0
    assert "candidate_id" in result.candidates.columns


def test_load_candidates_full_csv():
    result = load_candidate_file("outputs/run_001/candidates_full.csv")
    assert len(result.candidates) > 0


def test_missing_optional_columns_are_tolerated(tmp_path):
    path = tmp_path / "minimal.tsv"
    pd.DataFrame([{"candidate_id": "x1", "isomeric_smiles": "CCO"}]).to_csv(path, sep="\t", index=False)
    result = load_candidate_file(path)
    assert "canonical_smiles" in result.candidates.columns
    assert result.warnings


def test_qc_files_load_and_skip_when_absent(tmp_path):
    qc = load_qc_files("outputs/run_001")
    assert qc["generation_summary"] is not None
    empty = load_qc_files(tmp_path)
    assert empty["generation_summary"] is None


def test_discover_run_folders_lists_only_loadable_immediate_children(tmp_path):
    outputs = tmp_path / "outputs"
    valid_tsv = outputs / "run_tsv"
    valid_json = outputs / "run_json"
    invalid = outputs / "not_a_run"
    nested = outputs / "group" / "nested_run"
    for folder in [valid_tsv, valid_json, invalid, nested]:
        folder.mkdir(parents=True)
    (valid_tsv / "candidates_full.tsv").write_text("candidate_id\tisomeric_smiles\n", encoding="utf-8")
    (valid_json / "candidates_metadata.json").write_text("[]", encoding="utf-8")
    (invalid / "notes.txt").write_text("not candidates", encoding="utf-8")
    (nested / "candidates_full.csv").write_text("candidate_id,isomeric_smiles\n", encoding="utf-8")

    assert discover_run_folders(outputs) == sorted(
        [valid_tsv.as_posix(), valid_json.as_posix()],
        key=str.casefold,
    )
