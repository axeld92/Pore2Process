import numpy as np

from .predict import carrier_flux, predict, score, states


def test_refined_exchange_preserves_parent_flux_and_marginals():
    flow = np.array([.4, .6])
    parent = np.array([[.3, .1], [.1, .5]])
    receiving = parent.T/flow[:, None]
    labels = np.array([[0, 1], [1, 2], [2, 2], [2, 3], [3, 0]])
    joint = carrier_flux(labels, flow, receiving, 2)
    np.testing.assert_allclose(joint.sum(axis=0), [.2, .2, .3, .3], rtol=1e-11)
    np.testing.assert_allclose(joint.sum(axis=1), [.2, .2, .3, .3], rtol=1e-11)
    np.testing.assert_allclose(joint.reshape(2, 2, 2, 2).sum(axis=(1, 3)), parent,
                               rtol=1e-13, atol=0.)


def test_memory_preserves_survivors_in_no_contact_route():
    flow = np.array([1.])
    receiving = np.array([[1.]])
    inert = np.arange(1., 401.)
    times = np.column_stack((inert, inert))
    labels = states(times, times, flow, 2)
    first_survives = inert <= 200
    resolved = predict(labels, first_survives, flow, receiving, 2)
    parent = predict(labels//2, first_survives, flow, receiving, 1)
    result = score(resolved, labels, np.zeros(400, dtype=int), first_survives,
                   ~first_survives, flow)
    assert parent["survivor_conditional_capture"] == .5
    np.testing.assert_allclose(resolved["survivor_conditional_capture"], .25/200.5)
    assert result["passed"]
