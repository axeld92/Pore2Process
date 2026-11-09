"""Actual Julia geometry/translation contract and saved-data rejection."""
import json

import pytest
from pore2process_demo.bridge import runtime
from graded_support.geometry import audit
from graded_support.model import read_protocol,profile,load_case
from graded_support.run import requests,native_case

def test_actual_julia_ratio_translation_and_fixed_interface():
    p=read_protocol(); jl=runtime()
    for bulk in p["geometry"]["bulk_radii_m"]:
        x,r=profile(p,bulk)
        assert r[0]==r[1]==r[-2]==r[-1]==4e-6
        assert x[-1]==40e-6
        assert r[2]==r[3]==bulk
    base=load_case(jl,8e14,8e14,p,2.5e8)
    better=load_case(jl,4e14,8e14,p,2.5e8)
    P0=3.12/3600/1e5
    assert base["pressure_drop_Pa"]==pytest.approx(100/3.6e6*(1/P0+2.5e8))
    assert better["pressure_drop_Pa"]==pytest.approx(100/3.6e6*(.5/P0+2.5e8))
    assert base["non_support_impedance_Pa_s_m"]==better["non_support_impedance_Pa_s_m"]
    assert better["flow_m3_s"]==pytest.approx(100/3.6e6)
    assert better["shaft_power_W"]==pytest.approx(better["pressure_drop_Pa"]*better["flow_m3_s"]/.75)
    assert 1-better["shaft_power_W"]/base["shaft_power_W"]==pytest.approx(.5*base["support_fraction"])
    jobs,names=requests(p,1e-3,1000.)
    assert len(jobs)==24 and len({name for name,_ in jobs})==24
    for name in names:
        chosen=[inputs for case_id,inputs in jobs if case_id.startswith('validation/'+name+'/')]
        assert sorted(c['pressure_Pa'] for c in chosen)==[-10,5,10,10,20]

def test_native_geometry_disconnection_rejected(tmp_path):
    rows=tmp_path/'voxels.csv'
    rows.write_text('i,j,k,material\n0,0,0,3\n1,0,0,1\n2,0,0,4\n')
    raw={'dx_m':1.,'mesh_span_m':2.,'axial_intervals':2,'inlet_nodes':1,'outlet_nodes':1}
    assert audit(rows,raw)['six_neighbor_connected_components']==1
    rows.write_text(rows.read_text()+'1,1,1,1\n')
    with pytest.raises(RuntimeError,match='disconnected'):
        audit(rows,raw)

def test_native_cache_tamper_is_rejected_before_reuse(tmp_path,monkeypatch):
    origin=tmp_path/'old'/'case';origin.mkdir(parents=True)
    binary=tmp_path/'test-binary';binary.write_bytes(b'not executed')
    monkeypatch.setattr('graded_support.run.BINARY',binary)
    from coupling.run import sha
    from graded_support.run import MEASUREMENTS
    inputs={'resolution':24}
    record={'binary_sha256':sha(binary),'inputs':inputs,'artifact_sha256':{name:'0'*64 for name in MEASUREMENTS}}
    for name in MEASUREMENTS: (origin/name).write_text('changed')
    (origin/'measurement.json').write_text(json.dumps(record))
    with pytest.raises(RuntimeError,match='artifact changed'):
        native_case(tmp_path/'new','case',inputs,read_protocol(),tmp_path/'old')

def test_nested_dwsim_summary_is_hash_checked(tmp_path):
    from coupling.run import sha
    from graded_support.check import verify
    nested=tmp_path/'process'/'nominal'/'straight'/'summary.json'
    nested.parent.mkdir(parents=True);nested.write_text('{"value":1}\n')
    s={'status':'passed','protocol':{},'artifact_sha256':{str(nested.relative_to(tmp_path)):sha(nested)}}
    (tmp_path/'summary.json').write_text(json.dumps(s))
    # Equivalent JSON with different bytes must still fail its original hash.
    nested.write_text('{"value": 1}\n')
    with pytest.raises(RuntimeError,match='Artifact hash mismatch'):
        verify(tmp_path/'summary.json')

def test_profile_verification_catches_connected_wrong_shape(tmp_path):
    rows=tmp_path/'voxels.csv'
    lines=['i,j,k,material']+[f'{i},{j},{k},{3 if i==0 else 4 if i==2 else 1}' for i in range(3) for j in range(2) for k in range(2)]
    rows.write_text('\n'.join(lines)+'\n')
    raw={'dx_m':1.,'mesh_span_m':2.,'axial_intervals':2,'inlet_nodes':4,'outlet_nodes':4,
         'bulk_radius_m':1.,'lead_length_each_m':0.}
    p={'geometry':{'body_knots_fraction':[0,1],'body_length_m':2.,'neck_radius_m':1.,'radius_factors':[1,1]}}
    assert audit(rows,raw,p)['declared_profile_matches_material_labels'] is True
    # Removing an interior voxel preserves connectedness and both boundaries,
    # but no longer represents the prescribed circular sections.
    rows.write_text('\n'.join(line for line in lines if line!='1,1,1,1')+'\n')
    with pytest.raises(RuntimeError,match='declared radius profile'):
        audit(rows,raw,p)
