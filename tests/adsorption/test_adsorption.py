"""Independent unit, conservation, source and actual process contracts."""
from pathlib import Path
import numpy as np
import pytest
from examples.adsorption.bridge import column, module
from examples.adsorption.evidence import evidence
from examples.adsorption.check import artifacts, topology
from pore2process_demo.runtime.openlb import ROOT
from examples.adsorption.run import read
from examples.filtration.calibration.check import close_tree

def reference():return read(ROOT/'examples/adsorption/reference/adsorption_reference.json')

def test_source_units_and_unresolved_density_are_not_silently_fitted():
    e=evidence();d={x['name']:x for x in e['measurements_and_fitted_parameters']}
    assert d['Langmuir_affinity']['SI_value']==94.
    assert d['Langmuir_capacity']['SI_value']==pytest.approx(.04169)
    assert d['effective_intraparticle_diffusivity_at_50mg_L']['SI_value']==pytest.approx(6.5973e-12)
    assert e['source_consistency']['particle_density_if_mass_dimensions_porosity_are_used_kg_m3']!=pytest.approx(574.33)
    assert e['scenario_assumptions']['particle_envelope_dry_density_kg_m3']==1000.
    assert all(x is None for x in e['unavailable'].values())

def test_real_julia_column_replays_accepted_process_case():
    c=reference()['cases'][-1];r=column(**c['inputs'])
    close_tree({k:v for k,v in r.items() if not isinstance(v,list) and k!='total_s'},c['column_metrics'],'Real Julia adsorption replay')
    assert r['normalized_mass_residual']<1e-8
    assert np.asarray(r['final_loading_kg_kg']).shape==(640,)
    assert all(np.diff(r['x_m'])>0)
    assert r['adsorbent_mass_kg']==pytest.approx(4.)
    assert len(r['time_s'])==r['steps']+1

def test_measured_native_resistance_sets_actual_dwsim_load():
    c=reference()['cases'][-1];p=c['inputs'];a=c['dwsim'];i=a['inputs']
    load=module().pressure_load(c['permeability_m2'],a['inlet_viscosity_Pa_s'],p['length_m'],p['superficial_velocity_m_s'],p['area_m2'],.75)
    assert i['Outlet_pressure_Pa']-i['Inlet_pressure_Pa']==pytest.approx(float(load.pressure_drop_Pa),rel=1e-9)
    assert a['shaft_power_W']==pytest.approx(float(load.reference_shaft_power_W),rel=1e-8)
    base=Path(reference()['summary_path']).parent/f"dp{c['particle_diameter_m']*1e6:g}um/dwsim"
    assert len(topology(base/'water_pump.dwxmz'))==3

def test_particle_comparison_holds_inventory_and_flow_fixed():
    s=reference();c=s['cases']
    assert len({x['inputs']['superficial_velocity_m_s'] for x in c})==1
    assert len({x['column_metrics']['adsorbent_mass_kg'] for x in c})==1
    assert c[0]['inputs']['ldf_s_inv']/c[-1]['inputs']['ldf_s_inv']==pytest.approx(16.)
    assert s['verification']['experimental_column_validation_passed'] is False

def test_native_vtp_retains_si_fields_and_zero_based_ids():
    import pyvista as pv
    base=Path(reference()['summary_path']).parent
    mesh=pv.read(base/'native_bed_reference.vtp')
    assert mesh.n_points>0
    assert np.array_equal(mesh['point_id_zero_based'],np.arange(mesh.n_points))
    assert mesh['native_velocity_m_s'].shape==(mesh.n_points,3)
    assert mesh.field_data['coordinate_units'][0]=='m'

def test_unreached_breakthrough_is_json_null_compatible():
    r=column(length_m=1.,t_end_s=.01,dt_s=.01,cells=8)
    assert r['breakthrough_time_s'] is None and r['treated_volume_m3'] is None

def test_background_bridge_calls_are_rejected_before_entering_julia():
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=1) as pool:
        task=pool.submit(column,t_end_s=.01,dt_s=.01,cells=8)
        with pytest.raises(RuntimeError,match='main server thread'):task.result()

@pytest.mark.parametrize('kwargs',[{'porosity':1.},{'ldf_s_inv':-1.},{'cells':1},{'dt_s':float('nan')}])
def test_invalid_adsorption_input_fails_before_a_process_solve(kwargs):
    with pytest.raises(Exception,match='porosity|Porosity|ldf|cells|dt_s'):column(**kwargs)

def test_changed_saved_artifact_rejected(tmp_path):
    from pore2process_demo.io import sha
    f=tmp_path/'field.csv';f.write_text('C\n1\n');s={'artifacts_sha256':{'field.csv':sha(f)}}
    f.write_text('C\n2\n')
    with pytest.raises(ValueError,match='Changed artifact'):artifacts(tmp_path/'summary.json',s)
