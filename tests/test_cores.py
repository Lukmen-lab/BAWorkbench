from ba_conjugate_generator.cores import CORE_REGISTRY, get_cores


def test_core_registry_loads_all_seven():
    assert set(CORE_REGISTRY) == {"CA", "CDCA", "DCA", "LCA", "UDCA", "7oxoCA", "7oxoCDCA"}
    assert all(core.mol is not None for core in get_cores())


def test_each_core_has_expected_available_sites():
    assert CORE_REGISTRY["CA"].hydroxyl_sites == ("C3-OH", "C7-OH", "C12-OH")
    assert CORE_REGISTRY["CDCA"].hydroxyl_sites == ("C3-OH", "C7-OH")
    assert CORE_REGISTRY["DCA"].hydroxyl_sites == ("C3-OH", "C12-OH")
    assert CORE_REGISTRY["LCA"].hydroxyl_sites == ("C3-OH",)
    assert CORE_REGISTRY["UDCA"].hydroxyl_sites == ("C3-OH", "C7-OH")
    assert CORE_REGISTRY["7oxoCA"].hydroxyl_sites == ("C3-OH", "C12-OH")
    assert CORE_REGISTRY["7oxoCDCA"].hydroxyl_sites == ("C3-OH",)
