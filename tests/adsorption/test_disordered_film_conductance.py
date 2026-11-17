import numpy as np
import pytest

from examples.adsorption.disordered_film_conductance.model import (
    effective_coefficient, film_rate, observations, required_coefficient, spherical_area_by_half,
)
from examples.adsorption.disordered_film_conductance.run import transfer, wall_survival


def test_clean_film_recovers_exponential_removal():
    # For a constant removal rate h and plug residence T, the survival is
    # exp(-hT) and the alive occupation averages (1-exp(-hT))/h.
    n = 100_000
    h, volume = 0.37, 1.7e-6
    hit = -np.log1p(-(np.arange(n) + .5) / n) / h
    arrivals = np.tile(np.arange(1., 9.), (n, 1))
    hits = np.where((hit[:, None] <= arrivals) & (hit[:, None] > arrivals - 1), hit[:, None], np.nan)
    counts, occupation, _ = observations(np.column_stack((arrivals, hits)))
    volumes = np.array([volume / 2, volume / 2, volume])
    np.testing.assert_allclose(counts.sum(axis=0)[:2].sum(), counts.sum(axis=0)[-1])
    np.testing.assert_allclose(occupation.sum(axis=0)[:2].sum(), occupation.sum(axis=0)[-1])
    np.testing.assert_allclose(film_rate(counts.sum(axis=0), occupation.sum(axis=0), volumes),
                               h * volumes, rtol=5e-5)


def test_periodic_surface_caps_and_seam_clock():
    areas = spherical_area_by_half(np.array([0., 2.2]), 1., 4.)
    np.testing.assert_allclose(areas, np.pi * np.array([.8, 1.2]))
    arrivals = np.tile(np.arange(1., 9.), (2, 1))
    hits = np.full_like(arrivals, np.nan)
    hits[0, 4] = 4.
    counts, times, _ = observations(np.column_stack((arrivals, hits)))
    np.testing.assert_array_equal(counts, [[0, 1, 1], [0, 0, 0]])
    np.testing.assert_array_equal(times, [[4, 0, 4], [4, 4, 8]])
    hits[0, 4] = 3.
    with pytest.raises(ValueError, match="contact"):
        observations(np.column_stack((arrivals, hits)))


def test_scaling_keeps_series_resistance_and_mass_based_area():
    geo = dict(external_area_m2=[2.], voxel_void_m3=[5.])
    column = dict(external_area_m2=8., liquid_volume_m3=10., length_scale=1.,
                  velocity_scale=3., direct_fraction=.5, added_resistance_s_m=4000.)
    measured = 4e-5
    G = transfer(measured, geo, column, True)
    k = measured * column["velocity_scale"]
    assert G == pytest.approx(column["external_area_m2"] * effective_coefficient(k, .5, 4000.))
    assert required_coefficient(G, 8., .5, 4000.) == pytest.approx(k)
    assert transfer(measured, geo, column, False) == pytest.approx(
        column["external_area_m2"] * effective_coefficient(k / 2, .5, 4000.))


def test_wall_weighted_plateau_preserves_probability():
    stages = [
        (np.array([0, 0]), np.array([0, 1]), None, np.array([.4, .5])),
        (np.array([0, 1]), np.array([0, 0]), None, np.array([.9, .8])),
    ]
    assert wall_survival(stages, 3, 1) == pytest.approx((.4 * .9 + .5 * .8)**3 * .9)
