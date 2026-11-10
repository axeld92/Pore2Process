"""Independent geometry/physics contracts and saved real-simulator evidence."""
import copy
import json
from zipfile import ZipFile

import numpy as np
import pyvista as pv
import pytest

from pore2process_demo.bridge import runtime
from examples.filtration.star_channel.geometry import polygon, dimensions, write_stl
from examples.filtration.star_channel.model import translate, requests
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read

def evidence():
    r=read(ROOT/'examples/filtration/star_channel/reference/star_channel_reference.json')
    return r,r['native_for_replay'],r['TRT_support_for_replay'],r['protocol']

def test_declared_star_area_perimeter_and_stl(tmp_path):
    v=polygon();g=dimensions(v)
    assert v.shape==(12,2)
    assert g['open_area_m2']==pytest.approx(3*.0023*.0014)
    side=np.sqrt(.0023**2+.0014**2-2*.0023*.0014*np.cos(np.pi/6))
    assert g['wetted_perimeter_m']==pytest.approx(12*side)
    f=tmp_path/'star.stl';write_stl(f,v,-.001,.001,.0023)
    mesh=pv.read(f)
    assert mesh.n_open_edges==0
    assert mesh.volume==pytest.approx(g['open_area_m2']*.002,rel=1e-12)
    with pytest.raises(ValueError):polygon(inner=.003)

def test_julia_axial_profiles_and_conservative_biological_feed():
    _,n,s,p=evidence();jl=runtime()
    q=dict(support='straight',flux_LMH=12.,rejection=.65,viscosity_multiplier=1.,fixed_multiplier=1.)
    t=translate(jl,n,s,p,q);a=t['axial'];f=t['feed']
    assert all(np.asarray(a[k]).shape==(101,) for k in ('x_m','pressure_Pa','local_flux_m_s','channel_flow_m3_s'))
    assert a['mean_TMP_Pa']==pytest.approx(12/3.6e6*a['total_membrane_impedance_Pa_s_m'])
    assert a['hydraulic_diameter_m']==pytest.approx(4*n['geometry']['open_area_m2']/n['geometry']['wetted_perimeter_m'])
    assert a['inlet_Re']==20
    assert t['critical_flux_limit_LMH'] is None and t['local_flux_admissible'] is None
    assert np.all(np.diff(a['pressure_Pa'])<0) and np.all(np.diff(a['channel_flow_m3_s'])<0)
    assert np.all(np.diff(f['axial_COD_kg_m3'])>0)
    assert max(a['local_flux_m_s'])*3.6e6==pytest.approx(t['maximum_local_flux_LMH'])
    assert f['module_inlet_COD_kg_m3']==pytest.approx(.0687,rel=1e-14)
    assert 0<f['fresh_COD_kg_m3']<f['module_inlet_COD_kg_m3']<f['retentate_COD_kg_m3']
    assert 0<f['permeate_COD_kg_m3']<f['fresh_COD_kg_m3']
    assert max(abs(v) for k,v in f.items() if k.endswith('_residual_kg_s'))<1e-17
    b=translate(jl,n,s,p,{**q,'support':'widened_6um'})
    assert b['axial']['mean_TMP_Pa']<a['mean_TMP_Pa']
    assert b['process']['total_shaft_power_W']<t['process']['total_shaft_power_W']
    assert b['axial']['fixed_impedance_Pa_s_m']==a['fixed_impedance_Pa_s_m']

@pytest.mark.parametrize('change',[{'flux_LMH':100.},{'flux_LMH':float('nan')},{'flux_LMH':None},{'rejection':1.1},{'support':'unknown'},{'viscosity_multiplier':0.},{'inlet_Re':1388.}])
def test_invalid_inputs_and_range_fail_before_external_solver(change):
    _,n,s,p=evidence()
    q=dict(support='straight',flux_LMH=12.,rejection=.65,viscosity_multiplier=1.,fixed_multiplier=1.)
    with pytest.raises(ValueError):translate(runtime(),n,s,p,{**q,**change})

def test_saved_actual_dwsim_replays_new_trt_geometry():
    r,n,s,p=evidence();c=r['actual_case'];q=c['request']
    t=translate(runtime(),n,s,p,q)
    from examples.filtration.calibration.check import close_tree
    from examples.filtration.crossflow_recycle.model import validate
    close_tree(t,c['translation'])
    checked=validate(c['dwsim'],c['dwsim']['inputs'],t['process'],r['hydraulics'],s)
    assert all(checked['checks'].values())
    assert c['dwsim']['inputs']['Force_liquid'] is True
    bad=copy.deepcopy(c['dwsim'])
    bad['saturation_pressure_Pa']['Fresh']=bad['streams']['Fresh']['pressure_Pa']
    with pytest.raises(RuntimeError,match='saturation/freezing guard'):
        validate(bad,bad['inputs'],t['process'],r['hydraulics'],s)
    assert s['protocol']['numerics']['collision_model']=='TRT'
    assert all(c['raw']['collision_model']=='TRT' for c in s['native_cases'].values())
    assert n['cases']['grid48']['raw']['collision_model']=='TRT'

def test_corrupt_artifact_is_rejected(tmp_path):
    from pore2process_demo.io import sha
    from examples.filtration.star_channel.check import verify
    f=tmp_path/'state.json';f.write_text('{}')
    (tmp_path/'summary.json').write_text(json.dumps({'status':'passed','artifacts_sha256':{'state.json':sha(f)}}))
    f.write_text('{"changed":true}')
    with pytest.raises(ValueError,match='Changed process artifact'):verify(tmp_path/'summary.json')
