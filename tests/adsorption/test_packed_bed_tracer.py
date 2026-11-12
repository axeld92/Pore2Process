"""Analytical Julia checks and the curated native-to-macro data contract."""
from concurrent.futures import ThreadPoolExecutor
import json
import math
import numpy as np
import pytest
from pore2process_demo.paths import ROOT
from examples.adsorption.packed_bed_tracer.bridge import call
from examples.adsorption.column_contact.audit import frozen

def reference():
    return json.loads((ROOT/'examples/adsorption/packed_bed_tracer/reference/tracer_reference.json').read_text())

@pytest.mark.parametrize('Pe',[2.,10.,40.,100.])
def test_absorbing_boundary_moments_have_independent_formula(Pe):
    L=.008;v=.0005;D=v*L/Pe
    r=call('absorbing_moments',length_m=L,velocity_m_s=v,dispersion_m2_s=D)
    coth=1/math.tanh(Pe/2)
    assert r['mean_s']==pytest.approx(L/v*coth-2*D/v**2,rel=1e-10)
    variance=2*D*L/v**3*coth+(L/v/math.sinh(Pe/2))**2-8*D**2/v**4
    assert r['variance_s2']==pytest.approx(variance,rel=1e-10)
    p=call('derive_contact',length_m=L,mean_s=r['mean_s'],variance_s2=r['variance_s2'])
    assert p['effective_velocity_m_s']==pytest.approx(v,rel=1e-8)
    assert p['effective_dispersion_m2_s']==pytest.approx(D,rel=1e-8)

def test_exit_flux_contract_and_one_native_segment_boundary():
    s=reference();p=s['contact'];a=s['cases']['bed32'];m=a['moments']
    actual=call('derive_contact',length_m=a['transport_length_m'],mean_s=m['mean_s'],variance_s2=m['variance_s2'])
    for key,value in actual.items():assert p[key]==pytest.approx(value,rel=1e-10)
    assert p['calibration_case']=='bed32' and p['holdout_case']=='length32'
    assert p['inventory_substitution_applied'] is False
    assert p['effective_contact_volume_m3']==pytest.approx(a['raw']['flow_m3_s']*m['mean_s'])
    assert p['physical_fluid_volume_m3']==a['raw']['fluid_volume_m3']
    assert 'not experimental' in p['role']
    assert s['new_native_scalar_runs']==8 and s['new_DWSIM_runs']==0
    assert set(s['cases'])=={'uniform16','uniform32','bed16','bed24','bed32','time32','length32','diffusion24'}
    for c in s['cases'].values():
        assert c['raw']['collision']=='stock AdvectionDiffusionTRTdynamics'
        assert 'EXIT-FLUX' in c['scope']
        assert c['raw']['tau_even']==pytest.approx(.5+c['raw']['magic']/(c['raw']['tau_odd']-.5))
        expected=(c['raw']['tau_odd']-.5)*c['raw']['dx_m']**2/(3*c['raw']['diffusivity_m2_s'])
        assert c['raw']['dt_s']==pytest.approx(expected)
        assert c['normalized_mass_residual']<1e-8
    assert s['native_checks_passed']==all(all(c['checks'].values()) for c in s['cases'].values())

def test_ade_reference_has_conserved_physical_storage_and_copy_arrays():
    p=dict(length_m=.008,velocity_m_s=.0005,dispersion_m2_s=1e-7,t_end_s=40.,dt_s=.005,cells=400)
    r=call('solve_ade',**p)
    mass=np.array(r['inventory_per_area_m']);incoming=np.array(r['cumulative_inlet_per_area_m']);outgoing=np.array(r['cumulative_outlet_per_area_m'])
    assert np.max(abs(mass-incoming+outgoing))/incoming[-1]<1e-9
    assert np.array(r['final_concentration']).shape==(400,)
    assert np.array(r['final_concentration']).dtype==np.float64
    assert min(r['final_concentration'])>=0 and max(r['final_concentration'])<=1

def test_batch_component_and_prior_column_failures_remain_explicit():
    s=frozen()
    assert s['batch_fit']['diffusivity_m2_s']==1.846035489738998e-12
    assert s['physical_status']['batch_particle_size_comparison_passed']
    assert not s['physical_status']['full_transport_comparison_passed']

def test_julia_calls_stay_on_main_thread():
    call('absorbing_moments',length_m=1.,velocity_m_s=1.,dispersion_m2_s=.1)
    with ThreadPoolExecutor(max_workers=1) as pool:
        with pytest.raises(RuntimeError,match='main server thread'):
            pool.submit(call,'absorbing_moments',length_m=1.,velocity_m_s=1.,dispersion_m2_s=.1).result()

def test_unrepresentable_contact_shape_is_not_silently_fitted():
    with pytest.raises(Exception,match='representable'):
        call('derive_contact',length_m=1.,mean_s=1.,variance_s2=.5)
