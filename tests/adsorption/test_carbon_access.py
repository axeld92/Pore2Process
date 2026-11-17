import numpy as np
import pytest

from examples.adsorption.carbon_access.kernel import load
from examples.adsorption.disordered_memory_column.run import julia_models
from pore2process_demo.bridge import runtime


def test_isolation_and_equal_access_limits_and_signed_balance():
    augmented = load().solve_column
    _,standard = julia_models()
    jl = runtime()
    rules = jl.seval("[MemoryColumn.Delays([1,1],[1,2],[.4,.6],[.5,1.5])]")
    flow = np.array([.4,.6])
    contact = np.array([[.65],[.35]])
    p = dict(blocks=1,liquid_volume_m3=1.1e-6,dry_mass_kg=.001,flow_m3_s=1e-6,
             inlet_concentration_kg_m3=.15,t_end_s=50.,dt_s=.25,
             radial_cells=12,diffusivity_m2_s=2e-12,film_m_s=5e-5)
    full = standard(rules,flow,contact,**p)
    identical = augmented(rules,flow,contact,**p,direct_fraction=.4,slow_ratio=1.)
    isolated = augmented(rules,flow,contact,**p,direct_fraction=.4,slow_ratio=0.)
    slow = augmented(rules,flow,contact,**p,direct_fraction=.4,slow_ratio=.1)
    partial = standard(rules,flow,contact,**dict(p,dry_mass_kg=.4*p["dry_mass_kg"]))
    for name in ("adsorbed_inventory_kg","outlet_concentration_kg_m3"):
        np.testing.assert_allclose(np.array(getattr(identical,name)),np.array(getattr(full,name)),
                                   rtol=1e-9,atol=1e-14)
        np.testing.assert_allclose(np.array(getattr(isolated,name)),np.array(getattr(partial,name)),
                                   rtol=1e-9,atol=1e-14)
    assert np.all(np.array(isolated.adsorbed_inventory_by_cohort_kg)[1]==0.)
    assert np.array(slow.adsorbed_inventory_by_cohort_kg)[1,-1]>0.
    assert isolated.physical_adsorbent_mass_kg == slow.physical_adsorbent_mass_kg == p["dry_mass_kg"]
    assert isolated.equilibrium_adsorbed_ceiling_kg == pytest.approx(.4*slow.equilibrium_adsorbed_ceiling_kg)
    for r in (identical,isolated,slow):
        assert r.normalized_mass_residual < 1e-9
        assert r.normalized_interface_flux_mismatch < 1e-9
