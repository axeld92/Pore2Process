import numpy as np
import pytest

from examples.adsorption.disordered_early_bypass.model import (
    compose, plateau_bound, propagation, required_ntu, wall_stages,
)


def test_joint_contact_retains_timing_of_survivors():
    arrivals = np.array([[1.],[2.]])
    hits = np.array([[np.nan],[1.5]])
    stages = wall_stages(arrivals,hits,np.zeros_like(arrivals,dtype=int))
    pooled = propagation(stages["receiver"],1.)[0]
    joint = propagation(stages["joint"],1.)[0]
    np.testing.assert_allclose(pooled,[0.,.25,.25],atol=1e-14)
    np.testing.assert_allclose(joint,[0.,.5,0.],atol=1e-14)
    assert compose([joint],2,0).sum() == pytest.approx(.25)


def test_uniform_conductance_attains_lower_escape_bound():
    flow = np.array([.25,.75])
    n = 8
    total = 3.
    unequal = np.array([.8,.2])
    a = total/n*unequal/flow
    mixed_escape = (flow@(1/(1+a)))**n
    uniform = (flow@(1/(1+total/n)*np.ones(2)))**n
    assert mixed_escape > plateau_bound(total,n)
    assert uniform == pytest.approx(plateau_bound(total,n))
    assert plateau_bound(required_ntu(.05,n),n) == pytest.approx(.05)
