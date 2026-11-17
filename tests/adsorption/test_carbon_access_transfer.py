import numpy as np
import pytest

from examples.adsorption.carbon_access_transfer import kernel, source
from examples.adsorption.carbon_access_transfer.run import transferred_scenarios
from pore2process_demo.bridge import runtime


def test_series_resistance_stays_fixed_when_fluid_film_changes():
    base = transferred_scenarios(2e-5,2e-5)
    changed = transferred_scenarios(2e-5,4e-5)
    for name in base:
        assert base[name]["added_resistance_s_m"] == changed[name]["added_resistance_s_m"]
        assert base[name]["inputs"]["direct_fraction"] == changed[name]["inputs"]["direct_fraction"]
    assert changed["f50_slow10"]["inputs"]["slow_ratio"] == pytest.approx(1/19)
    assert changed["f50_isolated"]["inputs"]["slow_ratio"] == 0.


def test_heldout_operating_conditions_and_duplicate_trace():
    cases = source.operating_cases()
    data = source.markers()
    assert len(cases) == 6
    assert len(data) == 7
    a,b = (data[f"f{i}_q233_l25"] for i in (4,5))
    mids = [source.observed_crossing(x["time"],x["fraction"],.5) for x in (a,b)]
    assert abs(mids[0]-mids[1]) < a["dt"]+b["dt"]


def test_partial_period_preserves_water_holdup_and_outlet_routing():
    solver = kernel.solver()
    jl = runtime()
    rules = jl.seval("[MemoryColumn.Delays([1],[1],[1.],[.4]),MemoryColumn.Delays([1],[1],[1.],[.7])]")
    p = dict(blocks=1,terminal_sections=1,liquid_volume_m3=1.5e-6,dry_mass_kg=.001,flow_m3_s=1e-6,
             inlet_concentration_kg_m3=1.,t_end_s=3.,dt_s=.1,radial_cells=1,film_m_s=0.)
    out = solver(rules,np.array([1.]),np.ones((1,2)),**p)
    assert out.sections == 3
    assert out.delay_scale == pytest.approx(1.)
    assert out.physical_adsorbent_mass_kg == pytest.approx(.001)
    assert np.array(out.liquid_inventory_kg)[-1] == pytest.approx(1.5e-6)
    time = np.array(out.time_s)
    concentration = np.array(out.outlet_concentration_kg_m3)
    np.testing.assert_allclose(concentration[1:],time[1:]-.05>1.5,atol=1e-12)
    assert out.normalized_mass_residual<1e-12
