"""Controls change only the declared inputs; failed transport stays unavailable."""
import csv
import json

import pytest

from pore2process_demo.paths import ROOT
from examples.adsorption.transport_accuracy.diagnosis import zero_velocity_copy


def test_stationary_control_preserves_coordinates_density_and_node_ids(tmp_path):
    source = tmp_path/'source.csv'
    target = tmp_path/'stationary.csv'
    source.write_text('i,j,k,x_m,y_m,z_m,ux_m_s,uy_m_s,uz_m_s,density_kg_m3\n'
                      '0,1,2,1.25e-5,3e-5,4e-5,.01,-.002,.003,1000.012\n')
    before = source.read_bytes()
    zero_velocity_copy(source, target)
    assert source.read_bytes() == before
    with source.open() as stream:
        original = next(csv.DictReader(stream))
    with target.open() as stream:
        result = next(csv.DictReader(stream))
    for key in original:
        if key in ('ux_m_s', 'uy_m_s', 'uz_m_s'):
            assert float(result[key]) == 0.
        else:
            assert result[key] == original[key]


def test_actual_diagnosis_separates_diffusion_from_flow_and_matches_time_steps():
    p = ROOT/'examples/adsorption/transport_accuracy/reference/diagnosis_reference.json'
    s = json.loads(p.read_text())
    assert s['new_native_runs'] == 4
    assert not s['numerical_source_changed']
    assert s['cases']['no_advection']['raw']['flow_m3_s'] == 0.
    assert s['cases']['no_advection']['passed']
    a = s['cases']['half_dt']['raw']
    b = s['cases']['grid48_matched_dt']['raw']
    assert a['dt_s'] == pytest.approx(b['dt_s'], rel=1e-9)
    assert a['diffusivity_m2_s'] == b['diffusivity_m2_s'] == 1e-9
    assert not s['service_prediction_qualified']
    assert s['new_adsorption_predictions'] == s['new_DWSIM_runs'] == 0


def test_actual_native_python_handoff_audit_has_no_coordinate_or_velocity_mismatch():
    p = ROOT/'examples/adsorption/transport_accuracy/reference/diagnosis_reference.json'
    s = json.loads(p.read_text())
    for audit in s['field_audits'].values():
        assert audit['fluid_node_IDs_match']
        assert audit['physical_velocities_match']
        assert max(audit['coordinate_maximum_error_m'].values()) == 0.
        assert audit['native_maximum_relative_density_departure'] < 1e-4
        assert 'Not OpenLB population flux' in audit['proxy_definition']


def test_actual_relaxation_channels_recover_the_declared_physical_diffusivity():
    p = ROOT/'examples/adsorption/transport_accuracy/reference/diagnosis_reference.json'
    s = json.loads(p.read_text())
    for case in [s['original_case'], *s['cases'].values()]:
        r = case['raw']
        # Stock TRT omega2=1/(magic/(1/omega_even-.5)+.5).
        # The scalar diffusive first moment is in this odd channel.
        tau_odd_from_stock = .5+r['magic']/(r['tau_even']-.5)
        D = (tau_odd_from_stock-.5)*r['dx_m']**2/(3*r['dt_s'])
        assert tau_odd_from_stock == pytest.approx(r['tau_odd'], rel=1e-12)
        assert D == pytest.approx(1e-9, rel=1e-10)
