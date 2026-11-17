import numpy as np

from examples.adsorption.porosity_matched_film import handoff, native, run


def test_dense_packing_preserves_sphere_volume_and_dimensions():
    points, acceptance = native.packing()
    assert points.shape == (32, 3)
    assert acceptance > 0
    assert np.min(points) >= 0
    assert np.all(points < native.EXTENT)
    assert native.minimum_distance(points, native.EXTENT) > native.DIAMETER_M
    np.testing.assert_allclose(1 - native.SOLID_VOLUME / np.prod(native.EXTENT), .427, atol=1e-14)
    for n in (45, 60):
        np.testing.assert_allclose(np.array(native.grid_shape(n)) * native.SIDE / n, native.EXTENT)


def test_external_application_retains_stock_dynamics():
    original = (native.PACKING / "openlb/capillary3d.cpp").read_text()
    generated = native.application(original)
    for text in ("collision::TRT", "forcing::Kupershtokh", "IndicatorSphere3D<T>",
                 "boundary::BounceBack", "setPeriodicity({true,true,true})", "l.collideAndStream()"):
        assert generated.count(text) == original.count(text) > 0
    assert "c.nx()*h-h" in generated
    assert "if (c.n<16 || c.n%5!=0" in generated


def test_common_cut_uses_physical_location_at_each_grid():
    a, b = np.full(12, -1.), np.full(16, -1.)
    a[2], b[3] = 1., 2.
    assert handoff.find_common_cut([a, b]) == (1, 4)
    b[3] = np.inf
    assert handoff.find_common_cut([a, b]) is None


def test_refinement_monitors_both_escape_transfers():
    def case(G, area, liquid):
        return dict(clean_film=dict(carbon_area=dict(effective_conductance_m3_s=G,
            wall_weighted_escape=area), liquid_volume=dict(wall_weighted_escape=liquid)))

    coarse = case(3e-6, .035, .07)
    assert run.compare(coarse, case(3.06e-6, .04, .065))["passed"]
    assert not run.compare(coarse, case(3.2e-6, .04, .065))["passed"]
    assert not run.compare(coarse, case(3.06e-6, .04, .10))["passed"]
