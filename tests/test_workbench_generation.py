import pandas as pd

import ba_conjugate_workbench.generation_controls as generation_controls
from ba_conjugate_workbench.generation_controls import WorkbenchGenerationRequest, run_generation


def test_workbench_generation_includes_lca_and_tauro_conjugates(tmp_path):
    summary = run_generation(
        WorkbenchGenerationRequest(
            cores=["CA", "LCA"],
            mode="primary_docking",
            output=tmp_path,
        )
    )

    candidates = pd.read_csv(tmp_path / "candidates_full.tsv", sep="\t")
    csv_candidates = pd.read_csv(tmp_path / "candidates_full.csv")
    tauro = candidates[
        (candidates["partner_normalized_name"] == "taurine")
        & (candidates["reaction_channel"] == "C24_amide")
    ]
    unmodified = candidates[candidates["reaction_channel"] == "unmodified"]
    csv_unmodified = csv_candidates[csv_candidates["reaction_channel"] == "unmodified"]

    assert summary["candidate_count"] == len(candidates)
    assert summary["unmodified_ba_count"] == 2
    assert candidates["candidate_id"].is_unique
    assert set(candidates["parent_core_abbrev"]) == {"CA", "LCA"}
    assert set(tauro["parent_core_abbrev"]) == {"CA", "LCA"}
    assert set(tauro["molecular_formula"]) == {"C26H45NO7S", "C26H45NO5S"}
    assert tauro["site_validation_passed"].all()
    assert set(unmodified["parent_core_abbrev"]) == {"CA", "LCA"}
    assert set(unmodified["molecular_formula"]) == {"C24H40O5", "C24H40O3"}
    assert (unmodified["mass_shift_from_parent_core"] == 0.0).all()
    assert unmodified["docking_ready"].all()
    assert len(csv_unmodified) == 2
    assert set(csv_unmodified["parent_core_abbrev"]) == {"CA", "LCA"}
    assert csv_unmodified["candidate_name"].notna().all()
    assert csv_unmodified["isomeric_smiles"].notna().all()
    assert csv_unmodified["neutral_monoisotopic_mass"].notna().all()


def test_manual_review_filter_never_removes_unmodified_ba_references(tmp_path):
    summary = run_generation(
        WorkbenchGenerationRequest(
            cores=["CA", "LCA"],
            mode="manual_review",
            output=tmp_path,
            include_manual_review=False,
        )
    )

    candidates = pd.read_csv(tmp_path / "candidates_full.csv")
    unmodified = candidates[candidates["reaction_channel"] == "unmodified"]

    assert summary["candidate_count"] == 2
    assert summary["unmodified_ba_count"] == 2
    assert set(unmodified["parent_core_abbrev"]) == {"CA", "LCA"}
    assert not unmodified["manual_review_required"].any()
    assert unmodified["docking_ready"].all()


def test_workbench_repairs_missing_unmodified_rows_instead_of_crashing(tmp_path, monkeypatch):
    monkeypatch.setattr(
        generation_controls,
        "generate_candidates",
        lambda cores, partners, config: ([], [], []),
    )

    summary = run_generation(
        WorkbenchGenerationRequest(
            cores=["CA", "LCA"],
            mode="primary_docking",
            output=tmp_path,
        )
    )
    candidates = pd.read_csv(tmp_path / "candidates_full.csv")

    assert summary["candidate_count"] == 2
    assert summary["unmodified_ba_count"] == 2
    assert summary["unmodified_ba_repaired"]
    assert set(candidates["parent_core_abbrev"]) == {"CA", "LCA"}
    assert set(candidates["reaction_channel"]) == {"unmodified"}
