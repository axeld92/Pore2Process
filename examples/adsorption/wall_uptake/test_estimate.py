"""Check the uptake estimator against a constant-hazard passage."""
import numpy as np

from .estimate import uptake


def test_perfect_sink_balance_recovers_known_hazard():
    n = 100_000
    duration = np.full(n,.2)
    ntu = .7
    t_hit = -np.log1p(-(np.arange(n)+.5)/n)*duration/ntu
    sampled = np.where(t_hit<duration,t_hit,np.nan)
    result = uptake(np.zeros(n,dtype=int),duration,sampled,1)
    np.testing.assert_allclose(result["capture"],[1-np.exp(-ntu)],atol=1e-5)
    np.testing.assert_allclose(result["NTU"],[ntu],atol=1e-5)
    assert result["hit_count"][0]>0


def test_censoring_and_all_captured_have_finite_conductance():
    result = uptake(np.array([0,0,1,1]),np.ones(4),np.array([np.nan,np.nan,.1,.3]),2)
    np.testing.assert_allclose(result["NTU"],[0,5])
