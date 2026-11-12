import numpy as np

from .predict import observe, profiles


def test_reset_preserves_mass_and_only_changes_later_contact():
    flow = [.5, .5]
    joint = [[.5, 0.], [0., .5]]
    untouched = profiles([1., 1.], flow, joint, 2, periods=4, reset_every=2)
    repeat = profiles([1., 0.], flow, joint, 2, periods=4)
    reset = profiles([1., 0.], flow, joint, 2, periods=4, reset_every=2)
    np.testing.assert_allclose(untouched["survival"], [1., 1., 1., 1.])
    np.testing.assert_allclose(repeat["survival"], [.5, .5, .5, .5])
    np.testing.assert_allclose(reset["survival"], [.5, .5, .25, .25])


def test_open_trajectory_censor_bounds():
    nan = np.nan
    data = np.array([[1., 2., nan, nan], [1., 2., .5, nan], [1., nan, nan, nan],
                     [1., nan, nan, 1.2]])
    observed = observe(data)
    assert observed["recovery"] == [1., .5]
    assert observed["survival"] == [.75, .25]
    assert observed["survival_upper"] == [.75, .5]
    assert observed["maximum_unresolved_fraction"] == .25
