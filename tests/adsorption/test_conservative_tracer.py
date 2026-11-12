"""The selected handoff preserves measured Q and physical pore inventory."""
import json
import numpy as np
import pytest
from pore2process_demo.paths import ROOT
from examples.adsorption.packed_bed_tracer.conservative import derive
from examples.adsorption.packed_bed_tracer.bridge import call

def reference():
    return json.loads((ROOT/'examples/adsorption/packed_bed_tracer/reference/conservative_reference.json').read_text())

def test_only_tracer_dispersion_is_derived_with_no_inventory_or_velocity_fit():
    r=reference();d=r['derived']
    assert d['fitted_quantities']==['effective dispersion from short inert-tracer mean only']
    assert d['calibration_case']=='bed32' and d['holdout_case']=='length32'
    assert not d['velocity_fitted'] and not d['inventory_fitted'] and not d['adsorption_data_used']
    assert r['comparison_passed'] and r['numerical_refinement_passed']
    assert r['physical_inventory_substitution'] is False
    assert r['new_native_runs']==r['new_DWSIM_runs']==0

def test_actual_julia_rederives_D_from_short_native_data_only():
    n=json.loads((ROOT/'examples/adsorption/packed_bed_tracer/reference/tracer_reference.json').read_text())
    actual=derive(n['cases']['bed32']);saved=reference()['derived']
    assert actual==saved

@pytest.mark.parametrize('key',['bed32','length32'])
def test_julia_area_scaling_conserves_known_flow_liquid_volume_and_solute(key):
    ref=reference();r=ref['predictions'][key];p=r['inputs'];area=r['physical_fluid_area_m2']
    assert area*p['length_m']==pytest.approx(r['physical_fluid_volume_m3'],rel=1e-13)
    assert area*p['velocity_m_s']==pytest.approx(r['flow_m3_s'],rel=1e-13)
    assert p['dispersion_m2_s']==ref['derived']['effective_dispersion_m2_s']
    actual=call('solve_ade',**p)
    mass=area*np.array(actual['inventory_per_area_m']);feed=area*np.array(actual['cumulative_inlet_per_area_m']);exit=area*np.array(actual['cumulative_outlet_per_area_m'])
    assert np.max(abs(mass-feed+exit))/feed[-1]<1e-8
    metrics=r['metrics'];assert metrics['passed']
    assert metrics['steady_exit_flux_relative_error']<=.1
