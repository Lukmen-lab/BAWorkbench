# BA Conjugate Candidate Generator

Python/RDKit prototype for generating chemically defined bile-acid-derived candidate molecules for docking-guided mass spectrometry discovery.

The tool generates covalent bile acid candidates, exports isomeric SMILES and metadata for later docking, and reports exact-mass/formula/adduct data for LC-MS/MS follow-up. It does not run Boltz-2 or any docking pipeline.

## Scientific Purpose

This project supports hypothesis generation. It creates a broad but traceable library of bile acid conjugates, including known-style C24 amides and esters plus exploratory covalent products that may be useful for docking prioritization and untargeted MS/MS searching. Generated candidates should not be interpreted as evidence that a molecule is naturally present.

Use plausibility and warning labels such as `reference_known`, `literature_supported`, `microbial_plausible`, `chemically_plausible`, `exploratory`, `reactive_or_artifact_possible`, and `manual_review_required` to triage downstream work.

## Bile Acid Cores

The first prototype includes:

- `CA`: cholic acid, C3-OH, C7-OH, C12-OH, C24-COOH
- `CDCA`: chenodeoxycholic acid, C3-OH, C7-OH, C24-COOH
- `DCA`: deoxycholic acid, C3-OH, C12-OH, C24-COOH
- `LCA`: lithocholic acid, C3-OH, C24-COOH
- `UDCA`: ursodeoxycholic acid, C3-OH, C7-OH, C24-COOH
- `7oxoCA`: 7-oxo-cholic acid, C3-OH, C7-oxo, C12-OH, C24-COOH
- `7oxoCDCA`: 7-oxo-chenodeoxycholic acid, C3-OH, C7-oxo, C24-COOH

Core SMILES are validated by RDKit at load time. The LCA stereochemistry is checked against PubChem CID 9903; the 7-oxo structures are flagged with TODO notes for production-grade stereochemical review.

## Reaction Channels

Implemented covalent channels:

- `C24_amide`: BA-C24-COOH plus primary/secondary amines
- `C24_ester`: BA-C24-COOH plus alcohols or phenols
- `BA_O_acyl`: BA hydroxyl esterification with carboxylic acids/acyl-like partners
- `BA_O_glycoside`: BA hydroxyl glycosides from saccharide partners
- `BA_polymeric_fragment_model`: defined glucan/exopolysaccharide DP fragments
- `BA_O_sulfate`: covalent BA sulfate esters
- `BA_O_phosphate`: covalent BA phosphate esters
- `BA_reactive_carbonyl_adduct`: flagged hemiacetal/carbonyl adduct models

Excluded chemistry:

- BA sodium/potassium/ammonium salts
- magnesium, iron, manganese, calcium, or other metal complexes
- phosphate/sulfate salt complexes
- undefined inorganic adducts
- undefined polymers as exact molecules

Salt contexts are collapsed only when they provide a covalent organic partner. For example, sodium acetate contributes acetate/acetic acid chemistry, and calcium pantothenate contributes pantothenic acid chemistry.

Taurine is included as a built-in reference partner. The `C24_amide` channel therefore generates tauro conjugates for every selected core, including taurolithocholic acid for LCA.

## Install

Use an environment with RDKit available:

```bash
python -m pip install -e ".[test]"
```

If your platform cannot install RDKit from PyPI, use conda/mamba for RDKit and then install the package in editable mode.

## Run

Generate the broad default library:

```bash
python -m ba_conjugate_generator.cli generate \
  --cores CA CDCA DCA LCA UDCA 7oxoCA 7oxoCDCA \
  --output outputs/run_001 \
  --mode broad_covalent \
  --max-saccharide-dp 5
```

Summarize an output table:

```bash
python -m ba_conjugate_generator.cli summarize \
  --input outputs/run_001/candidates_full.tsv
```

Modes:

- `primary_docking`: C24 amides, C24 esters, small BA-O-acyl products, monosaccharide glycosides
- `broad_covalent`: all chemically defined covalent candidates, larger saccharides, vitamins, nucleobases, indole metabolites, reactive flagged candidates
- `manual_review`: speculative/reactive candidates are generated and flagged for review

## Outputs

Main exports:

- `candidates_full.tsv`
- `candidates_full.csv`
- `candidates_metadata.json`
- `smiles_for_docking.tsv`
- `boltz_ligands/` with simple per-ligand YAML records for later Boltz integration

QC exports:

- `invalid_candidates.tsv`
- `excluded_partners.tsv`
- `manual_review_candidates.tsv`
- `duplicate_or_isomeric_groups.tsv`
- `generation_summary.json`

Each generated candidate includes reaction channel, conjugation site, partner, parent core, molecular formula, neutral monoisotopic mass, expected adduct m/z values, docking readiness, plausibility labels, warning flags, and inclusion/exclusion rationale. Every selected BA core is also included once in its unmodified form (`candidate_category=unmodified_BA`, `reaction_channel=unmodified`) as a reference ligand.

Arginine C24 amides explicitly model and label attachment at the alpha-amino nitrogen and at both terminal guanidino nitrogens. The chain-connected guanidino nitrogen is not used. L-cysteine is included as the sulfur-containing amino-acid partner; disulfide-linked L-cystine is excluded from the default partner set. Symmetry-equivalent reaction paths are collapsed to one candidate row.

## Warning Flags

- `manual_review_required`: do not pass directly to docking without review
- `reactive_or_artifact_possible`: carbonyl/aldehyde or other reactive product model
- `large_flexible_warning`: large saccharides, vitamins, or flexible structures may be challenging for docking
- `polymeric_fragment_model`: a defined oligomer fragment represents a polymeric input
- `defined_degree_of_polymerization=DPn`: explicit glucan fragment size
- `anomeric_uncertainty`: glycoside stereochemistry is a model assumption
- `exclude_context_only`: input is an inorganic salt, metal complex, gas, or undefined context-only material

## Tests

```bash
pytest
```

The tests cover core loading, expected sites, representative amides/esters/acylations/glycosides, glucan fragments, salt exclusions, mass/adduct calculations, invalid structure reporting, exports, and CLI help.

## Site-Identity Validation And Migration

Bile-acid reaction sites are defined by persistent atom-map numbers in explicit mapped templates for all seven cores. Steroid hydroxyl channels resolve C3, C7, and C12 through those maps; C24 amide/ester chemistry resolves C24 independently. No BA position is assigned from RDKit atom order or substructure-match order.

Every generated candidate includes:

- `validated_actual_modified_position`
- `site_validation_passed`
- `site_validation_message`
- `structure_id`
- final-structure duplicate metadata

A candidate cannot be `docking_ready=True` when site validation fails.

Audit and repair an earlier export without overwriting it:

```bash
python -m ba_conjugate_generator.cli migrate-sites \
  --input outputs/run_001/candidates_full.csv \
  --output outputs/run_001_site_corrected
```

The migration retains the original candidate ID and SMILES, detects the original graph attachment site, moves the existing substituent to the requested mapped BA site, and writes corrected exports plus:

- `site_migration_audit.tsv`
- `site_migration_summary.json`

For the original 1,680-row run, the regression-checked migration finds 891 mismatches, repairs all 891, leaves zero post-repair mismatches, and causes no formula, exact-mass, or formal-charge changes.

## Bile Acid Conjugate Workbench

The Bile Acid Conjugate Workbench is a local Streamlit GUI layered on top of the BA Conjugate Candidate Generator. The generator remains the chemistry source of truth; the Workbench loads generator outputs, filters candidate tables, displays 2D/3D structures, supports curation, and exports docking/MS follow-up files. It does not implement docking.

Target platform: Windows 10 Enterprise, CPU-only, no GPU required.

### Workbench Setup

Create or update the mamba/conda environment:

```bat
setup_env.bat
```

The environment is named `ba_conjugate_workbench` and installs RDKit from conda-forge along with pandas, Streamlit, py3Dmol, PyYAML, OpenPyXL, PyArrow, Pillow, and pytest.

If `mamba` or `conda` is not on PATH, edit `setup_env.bat` and `run_workbench.bat` to call your Miniforge/Mambaforge/Anaconda `activate.bat`.

Launch the GUI:

```bat
run_workbench.bat
```

To stop the running server, click **Terminate Workbench** at the top of the
sidebar. The current Streamlit server shuts down cleanly, after which its
`mamba run` or `conda run` wrapper exits normally. The button does not stop
unrelated Python, mamba, or conda processes.

Manual launch:

```bat
mamba activate ba_conjugate_workbench
streamlit run ba_conjugate_workbench/app.py
```

### Loading Runs

Use the sidebar Project / Run controls to load an existing run folder such as:

```text
outputs/run_001
```

The Workbench looks for `candidates_full.tsv`, `candidates_full.csv`, or `candidates_metadata.json`. It also loads QC files when present:

- `invalid_candidates.tsv`
- `excluded_partners.tsv`
- `manual_review_candidates.tsv`
- `duplicate_or_isomeric_groups.tsv`
- `generation_summary.json`

Missing QC files are skipped gracefully.

### Filtering And Table Review

Sidebar filters include BA core, candidate category, partner, reaction channel, conjugation site, modified position, neutral mass, `[M-H]-`, `[M+H]+`, docking status, warning flags, formula contains N/S/P, manual-review status, and text search across candidate identifiers, names, partners, SMILES, and notes.

The table is designed for at least 10,000 rows. It keeps the workflow table-based and generates 3D conformers only for selected candidates.

### 2D And 3D Structure Inspection

Selecting a candidate updates the metadata and structure panels.

The 2D panel uses RDKit depictions and highlights, when metadata or heuristics allow:

- BA core
- partner-derived group
- modified/conjugated site
- newly formed bond

The 3D panel generates an RDKit conformer on demand, caches it in `cache/conformer_cache/`, and displays it in a 3Dmol.js/py3Dmol-style browser viewer.

Mouse controls:

- left mouse drag = rotate
- right mouse drag = translate / pan
- mouse wheel = zoom

Atom label modes:

- Off
- Element only
- Element + RDKit atom index
- BA conventional labels + atom indices, when metadata or heuristics are available

Use the sidebar Clear conformer cache button to remove cached SDF conformers.

### Curation

Curation is saved separately and does not modify the original `candidates_full.tsv`.

Statuses:

- `unreviewed`
- `include_for_docking`
- `exclude_from_docking`
- `manual_review`
- `high_priority_for_MS`
- `interesting_but_later`

Curation files are written under the selected curation folder, for example:

- `curation_state.tsv`
- `curated_candidates_full.tsv`
- `curated_smiles_for_docking.tsv`
- `manual_review_notes.tsv`

### Workbench Exports

Available exports include:

- `visible_candidates.tsv`
- `selected_candidates.tsv`
- `curated_candidates_full.tsv`
- `curated_smiles_for_docking.tsv`
- `curated_candidates_metadata.json`
- `selected_ligands.sdf`
- `manual_review_notes.tsv`
- `curation_state.tsv`

Docking SMILES format:

```text
candidate_id<TAB>isomeric_smiles<TAB>candidate_name
```

### Known Limitations

- RDKit 2D depictions are static; full mouse interaction is in the 3D panel.
- 3D conformers are computational approximations.
- Large/flexible molecules may require careful manual review.
- A single generated conformer may not represent the relevant docking conformation.
- Generated candidates are hypotheses for prioritization and MS follow-up, not evidence of biological presence.
