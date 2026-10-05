import pandas as pd

from ba_conjugate_workbench.curation import load_curation, save_curation, update_curation
from ba_conjugate_workbench.exports import export_docking_smiles, export_table


def test_curation_state_save_load(tmp_path):
    df = update_curation(pd.DataFrame(), "cand1", "include_for_docking", "looks good", "tester")
    path = tmp_path / "curation_state.tsv"
    save_curation(df, path)
    loaded = load_curation(path)
    assert loaded.iloc[0]["curation_status"] == "include_for_docking"


def test_export_selected_docking_smiles(tmp_path):
    df = pd.DataFrame([{"candidate_id": "cand1", "isomeric_smiles": "CCO", "candidate_name": "ethanol"}])
    out = export_docking_smiles(df, tmp_path / "smiles.tsv")
    text = out.read_text(encoding="utf-8")
    assert "cand1" in text
    assert "CCO" in text


def test_exports_do_not_overwrite_original_without_explicit_path(tmp_path):
    original = tmp_path / "candidates_full.tsv"
    original.write_text("original", encoding="utf-8")
    export_table(pd.DataFrame([{"candidate_id": "cand1"}]), tmp_path / "visible_candidates.tsv")
    assert original.read_text(encoding="utf-8") == "original"
