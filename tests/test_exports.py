from pathlib import Path

from ba_conjugate_generator.config import GeneratorConfig
from ba_conjugate_generator.cores import get_cores
from ba_conjugate_generator.exports import export_outputs
from ba_conjugate_generator.generator import generate_candidates
from ba_conjugate_generator.partners import built_in_partners


def test_export_files_are_created(tmp_path: Path):
    partners = [p for p in built_in_partners() if p.normalized_name in {"glycine", "acetic acid"}]
    rows, invalid, excluded = generate_candidates(get_cores(["CA"]), partners, GeneratorConfig(mode="primary_docking"))
    export_outputs(tmp_path, rows, invalid, excluded)
    for name in [
        "candidates_full.tsv",
        "candidates_full.csv",
        "candidates_metadata.json",
        "smiles_for_docking.tsv",
        "invalid_candidates.tsv",
        "excluded_partners.tsv",
        "manual_review_candidates.tsv",
        "duplicate_or_isomeric_groups.tsv",
        "generation_summary.json",
    ]:
        assert (tmp_path / name).exists()
    assert (tmp_path / "boltz_ligands").is_dir()
