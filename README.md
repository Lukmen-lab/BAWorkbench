# BAWorkbench

[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/Lukmen-lab/BAWorkbench/releases/tag/v0.1.0)
[![Tests](https://github.com/Lukmen-lab/BAWorkbench/actions/workflows/tests.yml/badge.svg)](https://github.com/Lukmen-lab/BAWorkbench/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**BAWorkbench v0.1** is an open-source Python/RDKit and Streamlit application for generating, inspecting, curating, and exporting chemically defined bile-acid conjugate candidates. It is intended for hypothesis generation before molecular docking and LC-MS/MS follow-up. It does not perform docking and does not establish that a generated molecule exists biologically.

**Use at your own risk. Lukmen Lab and the contributors do not assume any responsibility or liability for generated structures, chemical interpretations, docking results, experimental decisions, data loss, or any direct or indirect consequences of using this software. Users are solely responsible for validating structures and results before research, publication, or other use.** Read [DISCLAIMER.md](DISCLAIMER.md).

## Capabilities

- Generates candidates from CA, CDCA, DCA, LCA, UDCA, 7-oxo-CA, and 7-oxo-CDCA.
- Includes every selected unmodified bile acid as a reference ligand.
- Models C24 amides and esters, BA hydroxyl acylation and glycosylation, sulfate/phosphate esters, defined glucan fragments, and flagged reactive-carbonyl adducts.
- Exports isomeric/canonical SMILES, formulae, exact masses, common adduct m/z values, reaction and attachment metadata, and docking-readiness flags.
- Provides an interactive local Workbench for filtering, 2D/3D review, curation, and export.
- Validates named BA modification sites using persistent atom maps.

## Scientific scope and cautions

BAWorkbench creates traceable chemical hypotheses. Candidate plausibility labels and warning flags must be reviewed before docking or experimental interpretation. Generated 3D conformers are computational starting points, not validated bioactive conformations. The 7-oxo structures retain explicit source-review warnings in the code. See [LIMITATIONS.md](LIMITATIONS.md) before using results in a manuscript.

## Requirements

- Python 3.10 or newer (Python 3.11 is used by `environment.yml`)
- Windows, macOS, or Linux
- RDKit and the dependencies declared in [`pyproject.toml`](pyproject.toml)
- Mamba/Conda is recommended because it provides consistent RDKit installation

No GPU is required.

## Installation

Clone the repository:

```bash
git clone https://github.com/Lukmen-lab/BAWorkbench.git
cd BAWorkbench
```

### Recommended: Mamba or Conda

Windows:

```bat
setup_env.bat
run_workbench.bat
```

macOS/Linux:

```bash
chmod +x setup_env.sh run_workbench.sh
./setup_env.sh
./run_workbench.sh
```

### Pip

If RDKit wheels are supported on your platform:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
baworkbench
```

The app opens locally in a browser, normally at `http://localhost:8501`. Use **Terminate Workbench** at the top of the sidebar to stop its Streamlit server cleanly.

## Quick start

The Workbench initially loads the small versioned run under `examples/example_run`. Generate a complete default library from the sidebar, or use the command line:

```bash
ba-conjugate-generator generate \
  --cores CA CDCA DCA LCA UDCA 7oxoCA 7oxoCDCA \
  --output outputs/run_001 \
  --mode broad_covalent \
  --max-saccharide-dp 5
```

Windows PowerShell uses a backtick instead of `\` for multiline commands, or the command can be entered on one line.

Summarize an export:

```bash
ba-conjugate-generator summarize --input outputs/run_001/candidates_full.tsv
```

Generation modes:

- `primary_docking`: smaller, higher-priority covalent candidate set.
- `broad_covalent`: all chemically defined channels, including larger and flagged candidates.
- `manual_review`: emphasizes speculative/reactive products and review flags.

## Chemistry represented

Implemented channels include:

- `C24_amide`
- `C24_ester`
- `BA_O_acyl`
- `BA_O_glycoside`
- `BA_polymeric_fragment_model`
- `BA_O_sulfate`
- `BA_O_phosphate`
- `BA_reactive_carbonyl_adduct`

Inorganic salt complexes, undefined metal complexes, gases, and undefined polymers are excluded as exact candidates. Salt-labelled organic inputs are collapsed to their covalent organic chemistry where appropriate.

Arginine C24 amides explicitly distinguish the alpha-amino nitrogen and the two terminal guanidino nitrogens. L-cysteine is included; disulfide-linked L-cystine is not in the default partner set. Symmetry-equivalent reaction paths are deduplicated.

## Extending BAWorkbench

See [docs/extending.md](docs/extending.md) for worked guidance on adding BA cores, adding modification compounds, and deciding whether a new reaction requires code. User-supplied compounds can use existing reaction channels without changing Python; new BA cores and new reaction algorithms require mapped templates, implementation, and tests.

## Output files

Each generated run can contain:

- `candidates_full.tsv` and `candidates_full.csv`
- `candidates_metadata.json`
- `smiles_for_docking.tsv`
- `boltz_ligands/`
- `generation_summary.json`
- `invalid_candidates.tsv`
- `excluded_partners.tsv`
- `manual_review_candidates.tsv`
- `duplicate_or_isomeric_groups.tsv`

Runtime `outputs/` and `cache/` directories are intentionally excluded from Git.

## Workbench workflow

1. Select an existing run or configure and generate a new one.
2. Filter by BA core, category, partner, reaction channel, site, mass, docking status, or warning flags.
3. Select a row to inspect metadata and highlighted 2D/3D structures.
4. Assign a curation status and notes.
5. Export visible, selected, docking-ready, or manually reviewed candidates.

Curation is stored separately and does not modify the source candidate table.

## Reproducibility and tests

Run the test suite from the repository root:

```bash
python -m pytest -q
```

The suite covers structure generation, site identity, masses/adducts, deduplication, exports, migration, the Workbench loader and filters, UI interactions, and graceful shutdown. GitHub Actions runs the same suite for each push and pull request.

To record a manuscript analysis, retain the Git commit or release tag, generation mode, selected BA cores, maximum saccharide DP, input partner file (if any), and the generated `generation_summary.json`.

## Suggested methods text

> Candidate bile-acid conjugates were enumerated with BAWorkbench v0.1 (Lukmen Lab; release v0.1.0). Molecular structures were generated and validated with RDKit using explicitly mapped bile-acid reaction sites. Candidate tables included isomeric SMILES, molecular formulae, monoisotopic masses, common adduct m/z values, reaction metadata, and review flags. Generated structures were treated as hypotheses for docking and mass-spectrometric prioritization rather than evidence of biological occurrence.

Adapt this wording to the exact generation mode, release commit, filters, and downstream workflow used.

## Citation

Until a paper or archived DOI is available, cite the software repository and release:

> Lukmen Lab. *BAWorkbench* (version 0.1.0). GitHub. https://github.com/Lukmen-lab/BAWorkbench

GitHub and citation managers can read [`CITATION.cff`](CITATION.cff). Add the associated paper/DOI there after publication.

## Contributing and support

Bug reports and feature requests are welcome through [GitHub Issues](https://github.com/Lukmen-lab/BAWorkbench/issues). See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidance and [SECURITY.md](SECURITY.md) for private vulnerability reporting guidance.

## License

BAWorkbench is released under the [MIT License](LICENSE).
