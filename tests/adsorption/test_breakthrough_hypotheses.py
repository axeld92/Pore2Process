from pathlib import Path

import numpy as np
import pytest

from examples.adsorption.breakthrough_hypotheses.run import exchange_per_section, exposure_groups, solver
from pore2process_demo.bridge import runtime


def test_crossovers_retain_the_measured_relaxation_at_new_section_length():
    section = exchange_per_section([[.60, .40], [.10, .90]], .20, .004, .001)
    assert section["directed_crossover_flow_fractions"] == pytest.approx([.08, .08])
    assert section["projected_release_over_sample"] == pytest.approx(.50)
    assert (1 - section["mixing"])**4 == pytest.approx(.50)


def test_exposure_partition_uses_contact_per_unit_flow():
    flow = np.array([.1, .2, .3, .4])
    contact = np.array([[.15, .25], [.15, .25], [.55, .35], [.15, .15]])
    assert exposure_groups(flow, contact).tolist() == [2, 2, 2, 1]


def test_symmetric_regions_reproduce_the_same_column_and_conserve_phenol():
    model_path = Path(__file__).resolve().parents[2] / "examples/adsorption/breakthrough_hypotheses/model.jl"
    solve = solver(model_path)
    jl = runtime()
    rule = jl.MemoryColumn.Delays(np.array([1, 1]), np.array([1, 2]),
                                 np.array([.25, .75]), np.array([.1, .3]))
    p = dict(blocks=3, terminal_sections=0, flow_m3_s=1e-6, liquid_volume_m3=1e-6,
             dry_mass_kg=.002, inlet_concentration_kg_m3=.15, particle_diameter_m=.0008,
             particle_density_kg_m3=1200., diffusivity_m2_s=2e-12, film_m_s=1e-4,
             radial_cells=12, t_end_s=5., dt_s=.05, direct_fraction=.5, slow_ratio=.1)
    inputs = (jl.seval("rule -> [rule]")(rule), np.array([.25, .75]), np.array([[.4], [.6]]))
    control = solve(*inputs, *[np.array([1.])]*3, np.array([1, 1]), **p)
    split = solve(*inputs, *[np.array([.3, .7])]*3, np.array([1, 1, 2, 2]), mixing=.4, **p)
    for field in ("outlet_concentration_kg_m3", "liquid_inventory_kg", "adsorbed_inventory_kg"):
        np.testing.assert_allclose(np.array(getattr(split, field)), np.array(getattr(control, field)),
                                   rtol=1e-11, atol=1e-14)
    assert split.normalized_mass_residual < 1e-9
    assert split.physical_adsorbent_mass_kg == p["dry_mass_kg"]
