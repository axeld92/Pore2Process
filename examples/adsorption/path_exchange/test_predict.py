"""Check the equal-flow remapping and a limiting exchange case."""
import numpy as np

from .predict import equiprobable_delays, pair_prediction


def test_equal_flow_bins_preserve_mean():
    delays = np.array([.1, .3])
    weights = np.array([.2, .8])
    np.testing.assert_allclose(equiprobable_delays(delays, weights, 2), [.22, .3])


def test_rank_exchange_changes_spread_but_not_mean():
    n = 100
    delays = np.array([.1, .3])
    weights = np.array([.5, .5])
    times = np.linspace(0, 1, 10001)
    raw = np.linspace(.1, .3, n)
    _, persistent = pair_prediction(np.column_stack((raw, raw)), delays, weights, times)
    _, switched = pair_prediction(np.column_stack((raw, raw[::-1])), delays, weights, times)
    np.testing.assert_allclose([persistent["mean"], switched["mean"]], [.402, .402])
    np.testing.assert_allclose([persistent["variance"], switched["variance"]],
                               [.04 + .004**2 / 12, .004**2 / 12])
    np.testing.assert_allclose([persistent["fast_retention"], switched["slow_retention"]], [1., 0.])
