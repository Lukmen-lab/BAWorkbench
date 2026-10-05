from __future__ import annotations

import json

import pandas as pd
import pytest
from rdkit import Chem

from ba_conjugate_generator.cores import CORE_REGISTRY, get_cores
from ba_conjugate_generator.migration import migrate_previous_output
from ba_conjugate_generator.partners import built_in_partners
from ba_conjugate_generator.reactions import (
    ba_o_acyl,
    c24_amides,
    c24_esters,
    glycosides,
    reactive_carbonyl_adducts,
    resolve_steroid_sites,
    sulfate_phosphate_esters,
)
from ba_conjugate_generator.site_validation import clear_atom_maps, validate_modified_site


EXPECTED_SITES = {
    "CA": ("C3", "C7", "C12"),
    "CDCA": ("C3", "C7"),
    "DCA": ("C3", "C12"),
    "LCA": ("C3",),
    "UDCA": ("C3", "C7"),
    "7oxoCA": ("C3", "C12"),
    "7oxoCDCA": ("C3",),
}


def _partner(normalized_name):
    return next(p for p in built_in_partners() if p.normalized_name == normalized_name)


@pytest.mark.parametrize("core_key,positions", EXPECTED_SITES.items())
def test_mapped_core_templates_have_explicit_persistent_site_maps(core_key, positions):
    core = CORE_REGISTRY[core_key]
    mapped = core.mapped_mol
    maps = {atom.GetAtomMapNum() for atom in mapped.GetAtoms()}
    assert set(range(1, 25)).issubset(maps)
    assert {24, 124, 224}.issubset(maps)
    for position in positions:
        number = int(position[1:])
        assert number in maps
        assert 100 + number in maps


CHANNEL_CASES = [
    ("BA_O_acyl", lambda core, site: ba_o_acyl(core, _partner("acetic acid"), site)),
    ("BA_O_glycoside", lambda core, site: glycosides(core, _partner("glucose"), site)),
    ("BA_O_sulfation", lambda core, site: [p for p in sulfate_phosphate_esters(core, site) if p.reaction_channel == "BA_O_sulfation"]),
    ("BA_O_phosphorylation", lambda core, site: [p for p in sulfate_phosphate_esters(core, site) if p.reaction_channel == "BA_O_phosphorylation"]),
    ("reactive_carbonyl", lambda core, site: reactive_carbonyl_adducts(core, _partner("acrolein"), site)),
]


@pytest.mark.parametrize("channel,builder", CHANNEL_CASES)
@pytest.mark.parametrize("core_key,positions", EXPECTED_SITES.items())
def test_every_steroid_oh_channel_edits_the_explicit_requested_site(channel, builder, core_key, positions):
    core = CORE_REGISTRY[core_key]
    for position in positions:
        products = builder(core, position)
        assert products, f"No {channel} product for {core_key} {position}"
        for product in products:
            assert product.site.position == position
            validation = validate_modified_site(core, product.mol, position)
            assert validation.passed, validation.message
            clean = clear_atom_maps(product.mol)
            exported = Chem.MolToSmiles(clean, isomericSmiles=True, canonical=True)
            assert ":103" not in exported


@pytest.mark.parametrize("core_key", EXPECTED_SITES)
def test_c24_channels_only_modify_c24(core_key):
    core = CORE_REGISTRY[core_key]
    products = c24_amides(core, _partner("glycine")) + c24_esters(core, _partner("ethanol"))
    assert products
    for product in products:
        assert product.site.position == "C24"
        validation = validate_modified_site(core, product.mol, "C24")
        assert validation.passed, validation.message


def test_absent_requested_site_is_structured_skip_without_substitution():
    resolution = resolve_steroid_sites(CORE_REGISTRY["CDCA"], "C12")
    assert resolution.sites == ()
    assert "no steroid hydroxyl at C12" in resolution.skipped_reason
    assert ba_o_acyl(CORE_REGISTRY["CDCA"], _partner("acetic acid"), "C12") == []


def test_first_run_migration_acceptance_counts_and_invariants(tmp_path):
    output = tmp_path / "corrected"
    fixture = "tests/fixtures/legacy_candidates.csv"
    summary = migrate_previous_output(fixture, output)
    assert summary == {
        "input": fixture,
        "total_rows": 2,
        "mismatches_before_repair": 1,
        "corrected_rows": 1,
        "failed_rows": 0,
        "mismatches_after_repair": 0,
        "formula_changes": 0,
        "exact_mass_changes": 0,
        "formal_charge_changes": 0,
    }
    corrected = pd.read_csv(output / "candidates_full.csv")
    audit = pd.read_csv(output / "site_migration_audit.tsv", sep="\t")
    assert len(corrected) == len(audit) == 2
    assert corrected["site_validation_passed"].all()
    assert (audit["formula_unchanged"] == True).all()  # noqa: E712
    assert (audit["exact_mass_unchanged"] == True).all()  # noqa: E712
    assert (audit["formal_charge_unchanged"] == True).all()  # noqa: E712
    assert {"original_candidate_id", "original_isomeric_smiles", "detected_original_site", "correction_status"}.issubset(corrected.columns)
    assert {"structure_id", "final_structure_duplicate_count", "is_duplicate_final_structure"}.issubset(corrected.columns)
    disk_summary = json.loads((output / "site_migration_summary.json").read_text(encoding="utf-8"))
    assert disk_summary["mismatches_after_repair"] == 0
