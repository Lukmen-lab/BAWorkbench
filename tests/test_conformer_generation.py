from ba_conjugate_workbench.conformer_cache import cache_path_for, clear_cache, generate_conformer


def test_conformer_generation_for_small_candidate(tmp_path):
    result = generate_conformer("CCO", cache_dir=tmp_path)
    assert result.mol is not None
    assert result.mol_block
    assert result.cache_path and result.cache_path.exists()


def test_conformer_cache_save_load(tmp_path):
    first = generate_conformer("CCO", cache_dir=tmp_path)
    second = generate_conformer("CCO", cache_dir=tmp_path)
    assert first.cache_path == second.cache_path
    assert second.from_cache


def test_mmff_fallback_to_uff_message(tmp_path):
    result = generate_conformer("[Na+]", cache_dir=tmp_path)
    assert result.mol is not None or result.error
    assert result.error or any("UFF" in msg or "MMFF" in msg for msg in result.status_messages)


def test_invalid_smiles_handling(tmp_path):
    result = generate_conformer("not-a-smiles", cache_dir=tmp_path)
    assert result.mol is None
    assert result.error


def test_clear_cache(tmp_path):
    generate_conformer("CCO", cache_dir=tmp_path)
    clear_cache(tmp_path)
    assert not list(tmp_path.glob("*.sdf"))
