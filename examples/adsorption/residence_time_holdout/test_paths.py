"""Independent checks of pulse timing, moments and the nonnegative fit."""
import numpy as np

from .paths import moments, paths, pulse_cdf


def test_parallel_paths_with_a_square_pulse():
    times = np.array([0., 1., 1.1, 1.4, 2.1, 2.2, 3., 3.4, 4.])
    delays = np.array([1., 3.])
    weights = np.array([.75, .25])
    pulse = .4
    np.testing.assert_allclose(pulse_cdf(times, delays, weights, pulse),
                               [0., 0., .1875, .75, .75, .75, .75, 1., 1.], atol=1e-14)
    np.testing.assert_allclose(pulse_cdf(times, 2 * delays, weights, pulse),
                               [0., 0., 0., 0., .1875, .375, .75, .75, .75], atol=1e-14)
    np.testing.assert_allclose(moments(delays, weights, pulse), [1.7, .75 + pulse**2 / 12])
    np.testing.assert_allclose(moments(2 * delays, weights, pulse), [3.2, 3. + pulse**2 / 12])


def test_fit_of_a_known_path_distribution():
    dt, pulse = .001, .004
    times = dt * np.arange(801)
    delays = dt * np.arange(797)
    weights = (.8 * np.exp(-.5 * ((delays - .08) / .004)**2) / .004
               + .2 * np.exp(-.5 * ((delays - .18) / .025)**2) / .025)
    weights /= weights.sum()
    transit = float(np.dot(delays, weights))
    expected = pulse_cdf(times, delays, weights, pulse)
    fit_delays, fit_weights, _ = paths(expected, dt, pulse, transit)
    assert np.all(fit_weights >= 0)
    np.testing.assert_allclose([fit_weights.sum(), np.dot(fit_delays, fit_weights)],
                               [1., transit], atol=1e-13)
    np.testing.assert_allclose(pulse_cdf(times, 2 * fit_delays, fit_weights, pulse),
                               pulse_cdf(times, 2 * delays, weights, pulse), atol=2e-4)
