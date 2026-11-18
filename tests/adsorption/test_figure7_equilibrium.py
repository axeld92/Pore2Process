import pytest

from examples.adsorption.figure7_equilibrium.run import inputs, service


def test_equilibrium_transfer_preserves_the_column_and_resolves_time_units():
    original = dict(flow_m3_s=7.3e-7, dry_mass_kg=.007, liquid_volume_m3=4.06e-6,
                    film_m_s=4.54e-5, capacity_kg_kg=.1115, affinity_m3_kg=139.)
    profile = dict(capacity_kg_kg=.103, affinity_m3_kg=173.782, diffusivity_m2_s=1.846e-12)
    changed = inputs(dict(parameters=original), profile, (48, 4., .01))
    for key in ("flow_m3_s", "dry_mass_kg", "liquid_volume_m3", "film_m_s"):
        assert changed[key] == original[key]
    assert changed["capacity_kg_kg"] == .103
    assert original["capacity_kg_kg"] == .1115
    assert changed["startup_dt_s"] == .01 and changed["dt_s"] == 4.


def test_service_uses_seconds_and_distinguishes_unresolved_crossings():
    assert service(120., 7.3e-7, 3e-4) == pytest.approx(
        dict(time_s=120., volume_m3=8.76e-5, shaft_energy_J=.036))
    assert service(None, 7.3e-7, 3e-4) is None
    with pytest.raises(ValueError, match="positive flow"):
        service(120., 0., 3e-4)
