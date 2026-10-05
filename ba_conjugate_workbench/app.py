"""Bile Acid Conjugate Workbench Streamlit app."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from ba_conjugate_generator.cores import CORE_REGISTRY

from ba_conjugate_workbench.conformer_cache import cache_size, clear_cache, generate_conformer
from ba_conjugate_workbench.curation import CURATION_STATUSES, load_curation, merge_curation, save_curation, update_curation
from ba_conjugate_workbench.data_loader import discover_run_folders, load_candidate_file, load_qc_files, load_run_folder
from ba_conjugate_workbench.exports import export_docking_smiles, export_metadata_json, export_sdf, export_table
from ba_conjugate_workbench.filters import FilterSpec, apply_filters, column_options, warning_flag_options
from ba_conjugate_workbench.generation_controls import WorkbenchGenerationRequest, run_generation
from ba_conjugate_workbench.highlighting import large_flexible_warning
from ba_conjugate_workbench.process_control import schedule_workbench_shutdown
from ba_conjugate_workbench.qc_view import qc_summary
from ba_conjugate_workbench.state import init_state
from ba_conjugate_workbench.structure_2d import depiction_svg
from ba_conjugate_workbench.structure_3d import molblock_to_html
from ba_conjugate_workbench.table_view import DETAIL_COLUMNS, OPTIONAL_COLUMNS, candidate_id_from_selection, display_columns
from ba_conjugate_workbench.utils import DEFAULT_CACHE_DIR, candidate_smiles, ensure_directory


st.set_page_config(page_title="Bile Acid Conjugate Workbench", layout="wide")
init_state(st)


def _reset_candidate_selection() -> None:
    st.session_state.selected_candidate_id = ""
    st.session_state.last_detail_candidate_id = ""
    st.session_state.last_3d_result = None


def _run_folder_changed() -> None:
    st.session_state.candidate_file = ""
    _reset_candidate_selection()


def _candidate_file_changed() -> None:
    _reset_candidate_selection()


pending_run_folder = st.session_state.pop("pending_run_folder", "")
if pending_run_folder:
    st.session_state.run_folder = str(pending_run_folder)
    st.session_state.candidate_file = ""
    _reset_candidate_selection()


def _load_current() -> tuple[pd.DataFrame, list[str], Path | None]:
    try:
        if st.session_state.candidate_file:
            result = load_candidate_file(st.session_state.candidate_file)
        else:
            result = load_run_folder(st.session_state.run_folder)
        curation_path = Path(st.session_state.curation_folder) / "curation_state.tsv"
        curation = load_curation(curation_path)
        df = merge_curation(result.candidates, curation)
        st.session_state.curation_df = curation
        st.session_state.last_loaded_file = str(result.source_path)
        return df, result.warnings, result.source_path
    except Exception as exc:
        st.error(f"Could not load candidates: {exc}")
        return pd.DataFrame(), [], None


st.title("Bile Acid Conjugate Workbench")

with st.sidebar:
    if st.button(
        "Terminate Workbench",
        help="Stop this Workbench server and its mamba/conda run wrapper.",
        use_container_width=True,
    ):
        schedule_workbench_shutdown()
        st.success("Workbench is shutting down. You can close this browser tab.")
        st.stop()

    st.header("Project / Run")
    run_folder_options = discover_run_folders()
    current_run_folder = str(st.session_state.run_folder).strip()
    if current_run_folder and current_run_folder not in run_folder_options:
        run_folder_options.insert(0, current_run_folder)
    st.selectbox(
        "Load existing run folder",
        run_folder_options,
        key="run_folder",
        accept_new_options=True,
        filter_mode="fuzzy",
        placeholder="Select or enter a run folder",
        on_change=_run_folder_changed,
    )
    st.text_input("Load candidates_full.tsv/csv/json", key="candidate_file", on_change=_candidate_file_changed)
    st.text_input("Select output / curation folder", key="curation_folder")
    run_name = st.text_input("Run name", value="workbench_run")
    reload_clicked = st.button("Reload")
    if reload_clicked:
        _reset_candidate_selection()

    st.header("Generation")
    mode = st.selectbox("Generation mode", ["primary_docking", "broad_covalent", "manual_review"], index=1)
    core_label_to_key = {"7-oxo-CA": "7oxoCA", "7-oxo-CDCA": "7oxoCDCA", "CA": "CA", "CDCA": "CDCA", "DCA": "DCA", "LCA": "LCA", "UDCA": "UDCA"}
    selected_core_labels = st.multiselect("BA core selection", list(core_label_to_key), default=list(core_label_to_key))
    max_dp = st.slider("Max saccharide DP", min_value=1, max_value=10, value=5)
    include_reactive = st.checkbox("Include reactive candidates", value=True)
    include_large = st.checkbox("Include large flexible candidates", value=True)
    include_manual = st.checkbox("Include manual-review candidates", value=True)
    generation_output = st.text_input("Generation output folder", value=f"outputs/{run_name}")
    if st.button("Generate candidates"):
        with st.spinner("Running existing generator..."):
            req = WorkbenchGenerationRequest(
                cores=[core_label_to_key[label] for label in selected_core_labels],
                mode=mode,
                output=generation_output,
                max_saccharide_dp=max_dp,
                include_reactive=include_reactive,
                include_large_flexible=include_large,
                include_manual_review=include_manual,
            )
            summary = run_generation(req)
            st.success(
                f"Generated {summary['candidate_count']} candidates, including "
                f"{summary['unmodified_ba_count']} unmodified BAs, in {summary['output']}"
            )
            st.session_state.pending_run_folder = str(summary["output"])
            st.rerun()

df, load_warnings, source_path = _load_current()
for warning in load_warnings:
    st.warning(warning)

with st.sidebar:
    st.header("Filtering")
    spec = FilterSpec(
        cores=st.multiselect("Parent BA core", column_options(df, "parent_core_abbrev")),
        categories=st.multiselect("Candidate category", column_options(df, "candidate_category")),
        partners=st.multiselect("Partner molecule", column_options(df, "partner_name")),
        reaction_channels=st.multiselect("Reaction channel", column_options(df, "reaction_channel")),
        conjugation_sites=st.multiselect("Conjugation site", column_options(df, "conjugation_site")),
        modified_positions=st.multiselect("Modified position", column_options(df, "modified_position")),
        docking_statuses=st.multiselect("Docking status", column_options(df, "docking_status")),
        warning_flags=st.multiselect("Warning flags", warning_flag_options(df)),
        mass_range=tuple(st.slider("Neutral mass range", 0.0, 2500.0, (0.0, 2500.0))) if not df.empty else (None, None),
        mz_minus_h_range=tuple(st.slider("m/z range [M-H]-", 0.0, 2500.0, (0.0, 2500.0))) if not df.empty else (None, None),
        mz_plus_h_range=tuple(st.slider("m/z range [M+H]+", 0.0, 2500.0, (0.0, 2500.0))) if not df.empty else (None, None),
        formula_contains_n=st.checkbox("Formula contains N"),
        formula_contains_s=st.checkbox("Formula contains S"),
        formula_contains_p=st.checkbox("Formula contains P"),
        manual_review_required={"Any": None, "Yes": True, "No": False}[st.selectbox("Manual review required", ["Any", "Yes", "No"])],
        text_search=st.text_input("Text search"),
    )

    st.header("Visualization")
    show_2d = st.checkbox("Show 2D structure", value=True)
    show_3d = st.checkbox("Show 3D structure", value=True)
    auto_3d = st.checkbox("Auto-generate 3D on row selection", value=False)
    show_h = st.checkbox("Show hydrogens", value=True)
    label_mode = st.selectbox("Atom label mode", ["Off", "Element only", "Element + RDKit atom index", "BA conventional labels + atom indices"])
    st.checkbox("Highlight BA core", value=True)
    st.checkbox("Highlight partner group", value=True)
    st.checkbox("Highlight new bond", value=True)
    st.checkbox("Highlight modified site", value=True)
    if st.button("Clear conformer cache"):
        clear_cache(DEFAULT_CACHE_DIR)
        st.success("Conformer cache cleared.")
    count, size_bytes = cache_size(DEFAULT_CACHE_DIR)
    st.caption(f"Conformer cache: {count} files, {size_bytes / 1024:.1f} KB")

visible = apply_filters(df, spec) if not df.empty else df

top_cols = st.columns(6)
top_cols[0].metric("Total loaded", len(df))
top_cols[1].metric("Visible", len(visible))
top_cols[2].metric("Docking-ready", int((visible.get("docking_ready", pd.Series(dtype=bool)).astype(str).str.lower() == "true").sum()) if not visible.empty and "docking_ready" in visible else 0)
top_cols[3].metric("Manual review", int((visible.get("manual_review_required", pd.Series(dtype=bool)).astype(str).str.lower() == "true").sum()) if not visible.empty and "manual_review_required" in visible else 0)
mass_values = pd.to_numeric(visible.get("neutral_monoisotopic_mass", pd.Series(dtype=float)), errors="coerce") if not visible.empty else pd.Series(dtype=float)
top_cols[4].metric("Mass min", f"{mass_values.min():.2f}" if not mass_values.empty else "NA")
top_cols[5].metric("Mass max", f"{mass_values.max():.2f}" if not mass_values.empty else "NA")
st.caption(f"Current run folder: {st.session_state.run_folder} | Last loaded file: {st.session_state.last_loaded_file or 'none'}")

with st.expander("QC Summary", expanded=False):
    qc = load_qc_files(st.session_state.run_folder) if st.session_state.run_folder else {}
    st.json(qc_summary(qc, df))

st.subheader("Candidate Table")
optional = st.multiselect("Additional columns", [c for c in OPTIONAL_COLUMNS if c in visible.columns])
st.caption("Click a row to update Selected Candidate Details. Click a column header to sort.")
table_df = visible[display_columns(visible, optional)].copy().reset_index(drop=True) if not visible.empty else visible.copy()
table_event = st.dataframe(
    table_df,
    hide_index=True,
    width="stretch",
    height=520,
    key="candidate_table",
    on_select="rerun",
    selection_mode="single-row",
)
selected_candidate_id = candidate_id_from_selection(table_df, list(table_event.selection.rows))
visible_candidate_ids = set(visible["candidate_id"].astype(str)) if not visible.empty and "candidate_id" in visible else set()
if selected_candidate_id is not None:
    st.session_state.selected_candidate_id = selected_candidate_id
elif st.session_state.selected_candidate_id not in visible_candidate_ids:
    st.session_state.selected_candidate_id = str(visible.iloc[0]["candidate_id"]) if visible_candidate_ids else ""

selected_row = None
if st.session_state.selected_candidate_id and not df.empty:
    matches = df[df["candidate_id"].astype(str) == st.session_state.selected_candidate_id]
    if not matches.empty:
        selected_row = matches.iloc[0].to_dict()

st.subheader("Selected Candidate Details")
if selected_row is None:
    st.info("Select a candidate row to inspect structures and metadata.")
else:
    detail_cols = st.columns(2)
    details = {k: selected_row.get(k, "") for k in DETAIL_COLUMNS if k in selected_row}
    detail_cols[0].json(details)
    smiles = candidate_smiles(selected_row)
    if smiles:
        detail_cols[1].code(smiles, language="text")
    warning = large_flexible_warning(selected_row)
    if warning:
        st.warning(warning)

    left, right = st.columns(2)
    with left:
        st.markdown("**Highlighted 2D Structure**")
        if show_2d:
            svg, warnings = depiction_svg(selected_row, label_mode=label_mode)
            for warning_text in warnings:
                st.warning(warning_text)
            if svg:
                components.html(svg, height=440, scrolling=True)
                st.download_button("Download 2D SVG", svg, file_name=f"{selected_row['candidate_id']}.svg", mime="image/svg+xml")
    with right:
        st.markdown("**Interactive 3D Structure**")
        make_3d = st.button("Generate 3D conformer") or auto_3d
        if show_3d and make_3d:
            with st.spinner("Generating or loading 3D conformer..."):
                result = generate_conformer(smiles, DEFAULT_CACHE_DIR, show_hydrogens=show_h)
                st.session_state.last_3d_result = result
            if result.error:
                st.error(result.error)
            else:
                for msg in result.status_messages:
                    st.info(msg)
                components.html(molblock_to_html(result.mol_block, selected_row, label_mode=label_mode), height=560)
                st.caption("Mouse controls: left drag rotates, right drag translates/pans, mouse wheel zooms.")
        elif show_3d:
            st.info("Click Generate 3D conformer to create or load a cached conformer.")

    st.subheader("Curation")
    curation_folder = ensure_directory(st.session_state.curation_folder)
    curation_path = curation_folder / "curation_state.tsv"
    selected_detail_id = str(selected_row["candidate_id"])
    if st.session_state.last_detail_candidate_id != selected_detail_id:
        status = str(selected_row.get("curation_status", "unreviewed"))
        st.session_state.curation_status_input = status if status in CURATION_STATUSES else "unreviewed"
        st.session_state.reviewer_input = str(selected_row.get("reviewer", ""))
        st.session_state.curation_notes_input = str(selected_row.get("curation_notes", ""))
        st.session_state.last_detail_candidate_id = selected_detail_id
    cur_status = st.selectbox("Curation status", CURATION_STATUSES, key="curation_status_input")
    reviewer = st.text_input("Reviewer", key="reviewer_input")
    notes = st.text_area("Curation notes", height=120, key="curation_notes_input")
    if st.button("Save curation for selected"):
        curation_df = load_curation(curation_path)
        curation_df = update_curation(curation_df, str(selected_row["candidate_id"]), cur_status, notes, reviewer)
        save_curation(curation_df, curation_path)
        st.success(f"Curation saved to {curation_path}")

with st.sidebar:
    st.header("Export")
    out_dir = ensure_directory(st.session_state.curation_folder)
    selected_export = df[df["candidate_id"].astype(str) == st.session_state.selected_candidate_id] if st.session_state.selected_candidate_id and not df.empty else pd.DataFrame()
    if st.button("Export visible candidates"):
        export_table(visible, out_dir / "visible_candidates.tsv")
        st.success("Exported visible_candidates.tsv")
    if st.button("Export selected candidates"):
        export_table(selected_export, out_dir / "selected_candidates.tsv")
        st.success("Exported selected_candidates.tsv")
    if st.button("Export selected docking SMILES"):
        export_docking_smiles(selected_export, out_dir / "curated_smiles_for_docking.tsv")
        st.success("Exported curated_smiles_for_docking.tsv")
    if st.button("Export curation notes"):
        export_table(load_curation(out_dir / "curation_state.tsv"), out_dir / "curation_state.tsv")
        st.success("Exported curation_state.tsv")
    if st.button("Export manual review list"):
        manual = df[df.get("curation_status", "") == "manual_review"] if not df.empty and "curation_status" in df.columns else pd.DataFrame()
        export_table(manual, out_dir / "manual_review_notes.tsv")
        st.success("Exported manual_review_notes.tsv")
    if st.button("Export curated metadata JSON"):
        curated = df[df.get("curation_status", "") != "unreviewed"] if not df.empty and "curation_status" in df.columns else pd.DataFrame()
        export_metadata_json(curated, out_dir / "curated_candidates_metadata.json")
        st.success("Exported curated_candidates_metadata.json")
    if st.button("Export selected 3D SDF"):
        if len(selected_export) > 50:
            st.warning("Large SDF export requested. This may take a while.")
        path, messages = export_sdf(selected_export, out_dir / "selected_ligands.sdf", DEFAULT_CACHE_DIR)
        st.success(f"Exported {path}")
        for msg in messages:
            st.warning(msg)
