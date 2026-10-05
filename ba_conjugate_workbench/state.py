"""Streamlit session-state defaults."""

from __future__ import annotations


def init_state(st) -> None:
    defaults = {
        "run_folder": "outputs/run_001",
        "candidate_file": "",
        "curation_folder": "outputs/workbench_curation",
        "selected_candidate_id": "",
        "curation_df": None,
        "last_loaded_file": "",
        "last_3d_result": None,
        "last_detail_candidate_id": "",
        "curation_status_input": "unreviewed",
        "reviewer_input": "",
        "curation_notes_input": "",
    }
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)
