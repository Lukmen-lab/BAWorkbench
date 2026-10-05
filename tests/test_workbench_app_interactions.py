from streamlit.testing.v1 import AppTest


def test_candidate_row_selection_updates_selected_details_state():
    app = AppTest.from_file("ba_conjugate_workbench/app.py", default_timeout=60).run()
    assert not app.exception

    candidate_table = app.dataframe[0].value
    first_candidate = str(candidate_table.iloc[0]["candidate_id"])
    second_candidate = str(candidate_table.iloc[1]["candidate_id"])
    assert app.session_state["selected_candidate_id"] == first_candidate

    app.session_state["candidate_table"] = {
        "selection": {"rows": [1], "columns": []}
    }
    app.run(timeout=60)

    assert not app.exception
    assert app.session_state["selected_candidate_id"] == second_candidate
    assert app.session_state["last_detail_candidate_id"] == second_candidate


def test_run_folder_control_exposes_discovered_runs():
    app = AppTest.from_file("ba_conjugate_workbench/app.py", default_timeout=60).run()
    assert not app.exception

    assert app.sidebar.button[0].label == "Terminate Workbench"

    run_folder = next(
        widget for widget in app.selectbox if widget.label == "Load existing run folder"
    )
    assert "outputs/run_001" in run_folder.options
    assert "outputs/run_001_site_corrected" in run_folder.options
