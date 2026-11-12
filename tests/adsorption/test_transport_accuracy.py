"""An engineering diagnostic pass never supplies unqualified particle inputs."""
import json
import numpy as np
import pytest

from examples.adsorption.transport_accuracy.analysis import service_gate, refinement_metrics, remove_input_moments
from examples.adsorption.packed_bed_tracer.bridge import moments
from examples.adsorption.column_transport.bridge import qualified_dispersion
from pore2process_demo.paths import ROOT


def test_complete_screening_is_not_a_qualified_process_prediction():
    # Contract fixture only; no claim of an executed native trial.
    case = dict(engineering_checks={'kinetic_mass': True}, engineering_screen_passed=True)
    gate = service_gate({'screen': case}, ['screen'])
    assert gate['eligible_for_closure_assessment']
    assert gate['closure_qualified'] is False
    assert gate['service_prediction_qualified'] is False
    assert gate['service_time_s'] is None
    assert gate['energy_per_on_spec_volume_Wh_m3'] is None
    with pytest.raises(ValueError):
        qualified_dispersion(gate)


def test_missing_or_failed_screen_is_retained():
    bad = dict(engineering_checks={'concentration_nonnegative': False}, engineering_screen_passed=False)
    gate = service_gate({'bad': bad}, ['bad', 'missing_diameter'])
    assert gate['missing_cases'] == ['missing_diameter']
    assert gate['failed_native_screens'] == {'bad': ['concentration_nonnegative']}
    assert gate['eligible_for_closure_assessment'] is False


def test_exit_flux_moments_scale_with_hydraulic_time_not_signal_amplitude():
    # Independent exponential RTD: mean=tau and variance=tau^2.
    tau = 3.
    t = np.linspace(tau / 20000, 25 * tau, 40000)
    F = 1. - np.exp(-t / tau)
    a = moments(t, 7.3e-7 * F)
    b = moments(t, 7.3e-9 * F)
    assert a['mean_s'] == pytest.approx(tau, rel=1e-6)
    assert a['variance_s2'] == pytest.approx(tau ** 2, rel=1e-5)
    assert a['mean_s'] == pytest.approx(b['mean_s'], rel=1e-12)
    assert a['variance_s2'] == pytest.approx(b['variance_s2'], rel=1e-12)
    assert a['steady_exit_flux_kg_s'] / b['steady_exit_flux_kg_s'] == pytest.approx(100.)
    assert a['tail_second_moment_fraction'] < .001


def test_converged_mean_does_not_hide_unconverged_spreading():
    a = {'RTD': dict(mean_s=1., variance_s2=.1)}
    b = {'RTD': dict(mean_s=1.01, variance_s2=.2)}
    r = refinement_metrics(a, b)
    assert r['relative_changes']['mean_s'] < .05
    assert not r['passed']
    assert refinement_metrics(a, {'RTD': None}) == dict(available=False, passed=False)


def test_known_input_duration_is_removed_from_both_mean_and_variance():
    # A known exponential transport kernel convolved with a ramp derivative.
    # These are moments of the independent distributions, not a native fixture.
    tau = .2
    T = .1
    raw = dict(mean_s=tau+T/2, variance_s2=tau*tau+T*T/28,
               tail_area_fraction=.001, tail_second_moment_fraction=.002)
    corrected = remove_input_moments(raw, dict(mean_s=T/2, variance_s2=T*T/28))
    assert corrected['mean_s'] == pytest.approx(tau)
    assert corrected['variance_s2'] == pytest.approx(tau*tau)
    assert corrected['equivalent_sections'] == pytest.approx(1.)
    assert corrected['tail_area_fraction'] >= raw['tail_area_fraction']
    with pytest.raises(ValueError, match='nonpositive'):
        remove_input_moments(raw, dict(mean_s=1., variance_s2=1.))


def test_actual_native_evidence_separates_small_uniform_drift_from_front_failure():
    folder = ROOT / 'examples/adsorption/transport_accuracy/reference'
    s = json.loads((folder / 'screen_reference.json').read_text())
    assert s['actual_native_runs'] == 4
    for diameter in ('d800', 'd450'):
        constant = s['cases'][diameter + '_constant_base']
        front = s['cases'][diameter + '_front_base']
        assert constant['engineering_screen_passed']
        assert 1e-6 < constant['uniform_maximum_departure'] < .02
        assert front['minimum_concentration'] < -.3
        assert not front['engineering_screen_passed']
    assert s['decision']['service_prediction_qualified'] is False
    p = json.loads((folder / 'profile_reference.json').read_text())
    assert p['actual_native_runs'] == 4
    for diameter in ('d800', 'd450'):
        assert p['cases'][diameter + '_free']['passed']
        assert not p['cases'][diameter + '_porous']['passed']
    assert p['service_prediction_qualified'] is False
