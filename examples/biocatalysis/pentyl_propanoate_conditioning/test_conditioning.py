import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_conditioning.run import (
    INPUT, bead_storage, bed_mean_remaining, exchange_time, normalized_moments,
    reuse_proxies, scaling_cases, trace_stability,
)
from pore2process_demo.io import read_json


def test_finite_batch_storage_closes_water_inventory():
    for total in (0., .01, .42, 1., 10., 80.):
        free, q = bead_storage(total, 25.4, .025, .42)
        assert free + 25.4*q == pytest.approx(total, abs=1e-14)
        assert 0 <= free <= total
        assert 0 <= q < .025
        assert q == pytest.approx(.025*free/(.42+free))
    assert bead_storage(1., 25.4, 0., .42) == pytest.approx((1., 0.))
    with pytest.raises(ValueError):
        bead_storage(-1., 25.4, .025, .42)


def test_exchange_moments_remove_shared_pulse_and_recover_known_return_time():
    # A synthetic mobile clock with a known mean and variance, plus a
    # Poisson number of exponentially distributed delays at fixed T.
    mobile_mean, mobile_var, beta, tau = 16., 24., .35, 1.2
    pulse_mean, pulse_var = 3.5, 49/12
    b = {"mean_min": mobile_mean + pulse_mean, "variance_min2": mobile_var + pulse_var}
    a = {"mean_min": (1+beta)*mobile_mean + pulse_mean,
         "variance_min2": (1+beta)**2*mobile_var + 2*beta*mobile_mean*tau + pulse_var}
    result = exchange_time(b, a, pulse_mean, pulse_var)
    assert result["immobile_to_mobile_storage_ratio"] == pytest.approx(beta)
    assert result["return_time_min"] == pytest.approx(tau)
    assert result["physical_nonnegative_return_time"] is True
    upper = exchange_time(b, a, 0., 49/4)
    assert upper["return_time_min"] >= tau


def test_normalized_tracer_moments_do_not_depend_on_amplitude():
    t = np.linspace(0., 2., 1001)
    c = np.minimum(t, 2-t)
    rows = lambda scale: [{"time_min": ti, "concentration_mg_ml": ci*scale} for ti, ci in zip(t, c)]
    first, second = (normalized_moments(rows(scale)) for scale in (1., 5.))
    assert first["mean_min"] == pytest.approx(1.)
    assert first["variance_min2"] == pytest.approx(1/6, abs=1e-6)
    assert second["mean_min"] == pytest.approx(first["mean_min"])
    assert second["variance_min2"] == pytest.approx(first["variance_min2"])
    assert second["area_mg_min_ml"] == pytest.approx(5*first["area_mg_min_ml"])


def test_startup_contact_age_average():
    # An exact linear decay for the first ten minutes. At t=5 min, a
    # five-minute displacement exposes different sections for 0–5 min.
    t, remaining = [0., 600., 1200.], [1., 0., 0.]
    assert bed_mean_remaining(t, remaining, 5., 5.) == pytest.approx(.75)
    assert bed_mean_remaining(t, remaining, 0., 5.) == pytest.approx(1.)
    assert bed_mean_remaining(t, remaining, 20., 5.) == pytest.approx(0.)


def test_source_proxies_and_independent_illustrative_histories():
    inputs = read_json(INPUT)
    reuse = reuse_proxies(inputs)
    assert reuse["acetate_time_to_70pct_h"]["1"] == pytest.approx(4.54, abs=.05)
    assert reuse["acetate_time_to_70pct_h"]["10"] == pytest.approx(4.97, abs=.05)
    assert reuse["acetate_mild_time_scale"] == pytest.approx(.91, abs=.01)
    independent, illustration = scaling_cases(reuse, .64)
    assert independent["neat_acid_vs_alcohol_proxy"] == pytest.approx(1.7/46.2)
    assert independent["acid_1_5_reuse_proxy"] == pytest.approx(18.5/54.6)
    assert illustration["continued_decay_180min_illustration"] == pytest.approx(.64**6)
    assert illustration["recovered_wetting_180min_illustration"] == pytest.approx(1-.36**6)
    stability = trace_stability(inputs)
    assert stability["observed_30_to_180min_span_pp"] == pytest.approx(1.215, abs=.01)
