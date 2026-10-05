# Extending BAWorkbench

This guide describes the extension points in BAWorkbench v0.1. Test every extension with RDKit and review the scientific assumptions before using it in a paper or experiment.

## Add a modification compound without changing Python

The command-line generator accepts a JSON-compatible YAML partner file:

```bash
ba-conjugate-generator generate \
  --cores CA LCA \
  --partners my_partners.yaml \
  --mode primary_docking \
  --output outputs/my_run
```

The file is a list of partner records. Select one or more existing channels in `allowed_reaction_channels` and provide a valid isomeric SMILES:

```yaml
- partner_name: "Example alcohol"
  normalized_name: "example alcohol"
  synonyms: ["Example-OH"]
  smiles: "CCCO"
  formula: null
  source_category: "user_supplied"
  functional_group_tags: ["alcohol"]
  allowed_reaction_channels: ["C24_ester"]
  preferred_reaction_channels: ["C24_ester"]
  plausibility_class: "exploratory"
  warning_flags: ["manual_review_required"]
  notes: "Document the source and why this product is being screened."
  original_labels: ["Example alcohol"]
```

Use `C24_amide` for an amine, `C24_ester` for an alcohol or phenol, `BA_O_acyl` for a carboxylic-acid/acyl-like partner, `BA_O_glycoside` for a sugar partner, or `reactive_carbonyl` for an aldehyde/ketone model. The channel must match the functional-group assumptions enforced in `ba_conjugate_generator/reactions.py`; merely listing a channel does not make an invalid structure chemically valid.

The Workbench uses the built-in registry by default. For reproducible custom partner runs, retain the YAML file alongside the output and record its checksum or Git commit.

## Add a new BA core

A new BA core is a code change because the generator needs a mapped template and the Workbench needs a selectable label.

1. Add a `BileAcidCore` entry to `CORE_REGISTRY` in `ba_conjugate_generator/cores.py`.
2. Provide a validated unmodified `smiles` string and an explicit `mapped_smiles` template.
3. Preserve persistent atom maps: carbon positions use maps `1` `24`; hydroxyl oxygen at carbon `n` uses `100+n`; the C24 hydroxyl oxygen uses `124`; and the C24 carbonyl oxygen uses `224`.
4. Declare `sites` and `hydroxyl_sites` using the existing names such as `C3-OH`, `C7-OH`, `C12-OH`, and `C24-COOH`.
5. Add the display label and registry key to the Workbench selection mapping in `ba_conjugate_workbench/app.py`.
6. Add core-template, site-resolution, representative-reaction, formula/mass, and unmodified-reference tests.
7. Document the structure source, stereochemistry, and any uncertainty in the core notes and `LIMITATIONS.md`.

Example skeleton:

```python
"NEWBA": BileAcidCore(
    "NEWBA",
    "new bile acid",
    "<validated unmodified isomeric SMILES>",
    ("C3-OH", "C24-COOH"),
    ("C3-OH",),
    "Source: DOI or database identifier; stereochemistry reviewed on YYYY-MM-DD.",
    "<mapped template with persistent maps>",
),
```

Do not infer BA positions from RDKit atom order or a first substructure match. The site-validation tests are part of the contract.

## Add a reaction without changing Python

There is no general no-code reaction engine in v0.1. `data/reaction_rules.yaml` documents the existing channels; it is not a plug-in dispatcher. If a new compound fits an existing channel, add it as a partner as described above. This is the preferred no-code path.

## Add a genuinely new reaction channel

A new reaction algorithm requires code. At minimum:

1. Add a reaction-product builder to `ba_conjugate_generator/reactions.py`.
2. Add the channel to the generator dispatch in `ba_conjugate_generator/generator.py`.
3. Define candidate metadata, attachment-site semantics, bond type, warning flags, and validation behavior.
4. Add the channel description to `data/reaction_rules.yaml`, exports, and the README.
5. Add positive, invalid-structure, site-validation, deduplication, and mass/formula tests.
6. Confirm that unmodified BA reference rows remain present for every selected core.

Do not label a speculative product as `docking_ready` without documenting why. Keep exploratory or reactive products flagged for manual review.

## Before opening a pull request

Run:

```bash
python -m pytest -q
python -m pip wheel . --no-deps --no-build-isolation --wheel-dir .release-check
```

Include the exact input structures, sources, assumptions, expected reaction channel, and a short explanation of how the new tests demonstrate correctness.
