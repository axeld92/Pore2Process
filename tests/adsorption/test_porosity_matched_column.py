import numpy as np
import pytest

from examples.adsorption.porosity_matched_column.run import VARIANTS, curve_metrics, sensitivity


def test_grid_effects_hold_the_other_input_fixed():
    scores = {}
    for name, (film, transport) in VARIANTS.items():
        value = float(film - 2 * transport)
        scores[name] = dict(t05_s=dict(predicted_s=value), t50_s=dict(predicted_s=value),
            retained_kg_kg=value, common_dose_retained_kg_kg=value, rmse=value)
    result = sensitivity(scores)["t05_s"]
    assert result["paired_N45_minus_N60"] == 15.
    assert result["film_N45_minus_N60_fixed_transport"] == {"45": -15., "60": -15.}
    assert result["transport_N45_minus_N60_fixed_film"] == {"45": 30., "60": 30.}
    assert result["range_across_four_variants"] == [-75., -30.]


def test_crossing_and_retention_metrics_use_the_observation_window():
    basis = dict(dry_mass_kg=5., liquid_volume_m3=1., flow_m3_s=1.,
                 inlet_concentration_kg_m3=1., capacity_kg_kg=2., affinity_m3_kg=.1, feed_rate_kg_s=1.)
    saved = dict(inputs=basis, result=dict(
        time_s=[0., 10., 20.], outlet_concentration_kg_m3=[0., .5, 1.],
        liquid_inventory_kg=[0., 1., 1.], adsorbed_inventory_kg=[0., 4., 4.],
        adsorbed_inventory_by_cohort_kg=[[0., 3., 3.], [0., 1., 1.]],
        cumulative_outlet_kg=[0., 5., 15.],
        normalized_mass_residual=0., normalized_interface_flux_mismatch=0.))
    data = dict(time=np.array([10., 20.]), fraction=np.array([.5, 1.]))
    observed = dict(t05_s=1., t50_s=10., reading_t05_s=[.9, 1.1],
        reading_t50_s=[9., 11.], end_s=20., reading_retention_kg_kg=[.99, 1.01])
    metrics = curve_metrics(saved, basis, data, observed, common_time=5.)
    assert metrics["t05_s"]["predicted_s"] == pytest.approx(1.)
    assert metrics["t50_s"]["predicted_s"] == pytest.approx(10.)
    assert metrics["curve_screen_passed"] and metrics["retained_in_reading_range"]
    assert metrics["retained_kg_kg"] == pytest.approx(1.)
    assert metrics["common_dose_retained_kg_kg"] == pytest.approx(.5)
