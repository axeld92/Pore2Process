"""Actual Julia scale boundary and meaningful saved-evidence rejection checks."""
import copy
import json
import xml.etree.ElementTree as ET
import zipfile

import pytest
from pore2process_demo.bridge import runtime
from pore2process_demo.io import read_json as read; from pore2process_demo.paths import ROOT
from examples.filtration.crossflow_recycle.model import translate, adapter_inputs, topology, CONNECTIONS, validate
from examples.filtration.crossflow_recycle.run import jobs, check_protocol
from examples.filtration.crossflow_recycle.check import verify

def source():return read(ROOT/'examples/filtration/once_through/reference/process_filtration_reference.json')['source_for_Julia_replay']
def protocol():return read(ROOT/'examples/filtration/crossflow_recycle/protocol.json')

def test_real_julia_distinguishes_fresh_feed_module_and_return():
    s=source();p=protocol();m,c=translate(runtime(),s,s['supports'][0],p,.5,5.)
    qp=100/3.6e6;rho=s['fluid_anchor']['inlet_density_kg_m3']
    assert c['fresh_feed_reference_flow_m3_s']==pytest.approx(2*qp)
    assert c['module_reference_flow_m3_s']==pytest.approx(5*qp)
    assert c['recycle_reference_flow_m3_s']==pytest.approx(3*qp)
    assert c['purge_reference_flow_m3_s']==pytest.approx(qp)
    assert c['circulation_shaft_power_W']==pytest.approx((5*qp)**2*144e6/.75)
    assert c['feed_shaft_power_W']==pytest.approx(2*qp*m['pressure_drop_Pa']/.75)
    assert c['channel_pressure_drop_Pa']==pytest.approx(20000.)
    assert c['membrane_permeate_split_fraction']==pytest.approx(.2)
    assert c['purge_retentate_split_fraction']==pytest.approx(.25)
    request=adapter_inputs(p,s,c)
    assert request['Fresh_mass_flow_kg_s']==pytest.approx(rho*2*qp)
    assert request['Seed_module_mass_flow_kg_s']==pytest.approx(.5*rho*5*qp)

def test_scenario_matrix_keeps_recovery_and_circulation_independent():
    p=protocol();j=list(jobs(p,source()))
    assert len(j)==25 and sum(not v['validation_only'] for v in j)==21
    assert len({v['case_id'] for v in j})==25
    assert all(v['module_ratio']>1/v['recovery'] for v in j)
    assert sorted({v['module_ratio'] for v in j if v['condition'].startswith('circulation_')})==[10.,15.]
    assert {v['recovery'] for v in j if v['condition'].startswith('circulation_')}=={.5}

def test_changed_anchor_requires_review():
    p=protocol();p['temperature_K']=310.
    with pytest.raises(RuntimeError,match='changed anchor'):check_protocol(p,source())

def archive(path,edges):
    root=ET.Element('DWSIM_Simulation_Data');graphics=ET.SubElement(root,'GraphicObjects')
    sim=ET.SubElement(root,'SimulationObjects')
    for name in {v[0] for v in edges}|{v[2] for v in edges}:
        obj=ET.SubElement(graphics,'GraphicObject');ET.SubElement(obj,'Name').text='id-'+name;ET.SubElement(obj,'Tag').text=name
        outputs=[v for v in edges if v[0]==name and v[1]!='energy'];outs=ET.SubElement(obj,'OutputConnectors')
        for i in range(max((int(v[1]) for v in outputs),default=-1)+1):
            edge=next((v for v in outputs if int(v[1])==i),None)
            attrs={'IsAttached':'false'} if edge is None else {'IsAttached':'true','AttachedToObjID':'id-'+edge[2],'AttachedToConnIndex':edge[3]}
            ET.SubElement(outs,'Connector',attrs)
        for edge in edges:
            if edge[0]==name and edge[1]=='energy':ET.SubElement(ET.SubElement(obj,'EnergyConnector'),'Connector',{'IsAttached':'true','AttachedToObjID':'id-'+edge[2],'AttachedToConnIndex':edge[3]})
        if name in ('Tear-in','Circulation-in'):
            o=ET.SubElement(sim,'SimulationObject');ET.SubElement(o,'Name').text='id-'+name;ET.SubElement(o,'OverrideSingleCompoundFlashBehavior').text='true'
    with zipfile.ZipFile(path,'w') as z:z.writestr('saved.xml',ET.tostring(root))

def test_saved_recycle_and_heat_connector_contract(tmp_path):
    f=tmp_path/'saved.dwxmz';archive(f,CONNECTIONS);assert len(topology(f))==26
    archive(f,{v for v in CONNECTIONS if v[0]!='REC-1'})
    with pytest.raises(RuntimeError,match='topology differs'):topology(f)
    archive(f,{v for v in CONNECTIONS if v[1]!='energy'})
    with pytest.raises(RuntimeError,match='topology differs'):topology(f)

def test_changed_cooling_artifact_fails_before_replay(tmp_path):
    from pore2process_demo.io import sha
    f=tmp_path/'process/case/summary.json';f.parent.mkdir(parents=True);f.write_text('{"cooling_W":2}\n')
    s={'status':'passed','artifact_sha256':{'process/case/summary.json':sha(f)}}
    (tmp_path/'summary.json').write_text(json.dumps(s));f.write_text('{"cooling_W":0}\n')
    with pytest.raises(RuntimeError,match='Artifact hash mismatch'):verify(tmp_path/'summary.json')

def test_actual_saved_recycle_states_match_the_real_julia_model():
    r=read(ROOT/'examples/filtration/crossflow_recycle/reference/crossflow_recycle_reference.json');s=r['source_for_Julia_replay'];p=r['protocol']
    case=r['nominal_actual_DWSIM_case'];_,c=translate(runtime(),s,s['supports'][0],p,.5,5.)
    checks=validate(case['dwsim'],adapter_inputs(p,s,c),c,p,s)
    assert all(checks['checks'].values())
    assert checks['total_shaft_power_W']==pytest.approx(4.457454462461781)
    assert checks['circulation_shaft_fraction']==pytest.approx(.8308972567493434)
    assert checks['specific_shaft_energy_kWh_m3_reference_volume']==pytest.approx(.04457454462511821)

@pytest.mark.parametrize('change',['false_convergence','lost_mass','lost_heat','wrong_density'])
def test_passed_flags_cannot_hide_invalid_recycle_states(change):
    r=read(ROOT/'examples/filtration/crossflow_recycle/reference/crossflow_recycle_reference.json');case=r['nominal_actual_DWSIM_case'];a=copy.deepcopy(case['dwsim'])
    if change=='false_convergence':a['recycle']['converged']=False
    if change=='lost_mass':a['streams']['Return']['mass_flow_kg_s']*=.99
    if change=='lost_heat':
        a['powers']['cooling_removed_W']=a['powers']['cooling_energy_stream_W']=0.
    if change=='wrong_density':a['streams']['Circulation-in']['density_kg_m3']*=1.01
    with pytest.raises((RuntimeError,ValueError)):
        validate(a,a['inputs'],case['julia_process'],r['protocol'],r['source_for_Julia_replay'])
