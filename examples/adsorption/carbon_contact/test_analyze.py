"""Independent checks for fixed mass and shared wall accounting."""
import numpy as np

from .analyze import overlap, partition


def test_partitions_mass_once_and_preserves_unresolved_inventory():
    exposures = np.array([[2., 1.], [0., 4.], [0., 0.]])
    mass, unresolved = partition(exposures,np.array([3., 2., 5.]))
    np.testing.assert_allclose(mass,[2.,3.])
    assert unresolved==5.
    assert mass.sum()+unresolved==10.


def test_exclusive_area_bound_accounts_for_visits_from_other_classes():
    exposures = np.array([[1., 1.], [0., 1.], [1., 0.], [0., 0.]])
    result = overlap(exposures,np.ones(4))
    assert result["area_visited"]==.75
    assert result["shared_slow_and_other_area"]==.25
    assert result["maximum_potentially_slow_exclusive_area"]==.5
    assert result["slow_face_area_allocation_bounds"]==[.375,.625]
