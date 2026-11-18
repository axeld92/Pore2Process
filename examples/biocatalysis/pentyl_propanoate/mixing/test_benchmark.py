import numpy as np

from examples.biocatalysis.pentyl_propanoate.mixing.benchmark import moments, pulse_response


def test_rectangular_pulse_is_applied_in_physical_time():
    ages = np.full(32, 15.0)
    width = 7.0
    result = moments(ages, width)
    assert result["recovery"] == 1.0
    assert result["pulse_free_mean_min"] == 15.0
    assert result["pulse_free_variance_min2"] == 0.0
    assert result["rectangular_pulse_mean_min"] == 18.5
    assert result["rectangular_pulse_variance_min2"] == width**2 / 12
    e, f = pulse_response(np.array([14.0, 17.0, 20.0, 23.0]), ages, width)
    np.testing.assert_allclose(e, [0, 1 / width, 1 / width, 0])
    np.testing.assert_allclose(f, [0, 2 / width, 5 / width, 1])
