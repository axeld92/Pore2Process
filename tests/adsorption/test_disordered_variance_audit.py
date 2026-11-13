import numpy as np

from examples.adsorption.disordered_exchange.variance_audit import seam


def test_within_pair_covariance_isolated_at_fixed_exchange_and_delay():
    # Every class stays in place and has two ages with equally likely delays.
    labels = np.repeat(np.arange(10),2)
    age = np.tile([1.,3.],10)
    delay = age.copy()
    counts = 2*np.eye(10)
    result = seam(age,delay,labels,labels,counts,2*np.eye(10),5*np.eye(10),1.)
    np.testing.assert_allclose([result[k] for k in ("exchange","local_delay","propagation","memory_sort")],0.)
    assert result["memory_within_pair"] == result["memory"] == 2.
    assert result["measured_gain"] == 3.


def test_age_sorting_isolated_when_pair_delays_are_fixed():
    donor = np.repeat(np.arange(10),2)
    receiver = (donor+np.tile([0,1],10))%10
    age = np.tile([1.,3.],10)
    delay = age.copy()
    counts = np.bincount(10*donor+receiver,minlength=100).reshape(10,10)
    first = np.zeros((10,10))
    first[donor,receiver] = delay
    result = seam(age,delay,donor,receiver,counts,first,first**2,1.)
    np.testing.assert_allclose([result[k] for k in ("exchange","local_delay","propagation","memory_within_pair")],0.)
    assert result["memory_sort"] == result["memory"] == 2.
