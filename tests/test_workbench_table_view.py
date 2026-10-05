import pandas as pd

from ba_conjugate_workbench.table_view import candidate_id_from_selection


def test_candidate_id_from_native_single_row_selection():
    table = pd.DataFrame(
        [
            {"candidate_id": "first", "mass": 200.0},
            {"candidate_id": "second", "mass": 100.0},
        ]
    )

    assert candidate_id_from_selection(table, [1]) == "second"
    assert candidate_id_from_selection(table, []) is None
    assert candidate_id_from_selection(table, [10]) is None


def test_candidate_selection_follows_the_displayed_sorted_table():
    table = pd.DataFrame(
        [
            {"candidate_id": "first", "mass": 200.0},
            {"candidate_id": "second", "mass": 100.0},
        ]
    ).sort_values("mass").reset_index(drop=True)

    assert candidate_id_from_selection(table, [0]) == "second"
