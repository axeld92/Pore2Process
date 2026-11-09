"""Real Julia closure, flow contracts and rejection of altered saved evidence."""
import copy
import json
import xml.etree.ElementTree as ET
import zipfile

import pytest
from graded_support.model import read_protocol
from pore2process_demo.bridge import runtime
from process_filtration.model import translate,adapter_inputs,topology,CONNECTIONS,validate
from process_filtration.run import ROOT,read,check_protocol
from process_filtration.check import verify

def source_fixture():
    return {'protocol':read_protocol(),
        'literature_reference':{'fit':{'fixed_impedance_Pa_s_m':2.5e8}},
        'supports':[{'name':'straight','body_resistance_Pa_s_m3':8e14},
                    {'name':'widened','body_resistance_Pa_s_m3':4e14}],
        'fluid_anchor':{'inlet_density_kg_m3':1000.,'inlet_viscosity_Pa_s':1e-3}}

def test_actual_julia_full_feed_load_and_reference_mass_conversion():
    s=source_fixture();p=read(ROOT/'process_filtration/protocol.json');jl=runtime()
    check_protocol(p,s)
    a,b=translate(jl,s,s['supports'][1],.25,2e4)
    expected_tmp=100/3.6e6*(.5/(3.12/3600/1e5)+2.5e8)
    assert a['pressure_drop_Pa']==pytest.approx(expected_tmp)
    assert b['feed_reference_flow_m3_s']==pytest.approx(4*100/3.6e6)
    assert b['retentate_reference_flow_m3_s']==pytest.approx(3*100/3.6e6)
    assert b['shaft_power_W']==pytest.approx(4*100/3.6e6*(expected_tmp+2e4)/.75)
    assert b['specific_shaft_energy_kWh_m3']==pytest.approx((expected_tmp+2e4)/(.25*.75*3.6e6))
    request=adapter_inputs(p,s,b)
    assert request['Feed_mass_flow_kg_s']==pytest.approx(1000*4*100/3.6e6)
    assert request['TMP_Pa']==pytest.approx(expected_tmp)
    assert b['membrane_dissipation_W']+b['channel_dissipation_W']+b['retentate_throttle_dissipation_W']==pytest.approx(b['hydraulic_power_W'])

def test_changed_anchor_requires_new_native_property_review():
    s=source_fixture();p=read(ROOT/'process_filtration/protocol.json')
    s['protocol']['process']['temperature_K']=310.
    with pytest.raises(RuntimeError,match='changed anchor'): check_protocol(p,s)

def archive(path,connections):
    root=ET.Element('DWSIM_Simulation_Data');graphics=ET.SubElement(root,'GraphicObjects')
    names={v[0] for v in connections}|{v[2] for v in connections}
    for name in names:
        g=ET.SubElement(graphics,'GraphicObject');ET.SubElement(g,'Name').text='id-'+name;ET.SubElement(g,'Tag').text=name
        outs=ET.SubElement(g,'OutputConnectors')
        edges=sorted(v for v in connections if v[0]==name)
        if edges:
            for port in range(max(v[1] for v in edges)+1):
                e=next((v for v in edges if v[1]==port),None)
                ET.SubElement(outs,'Connector',{'IsAttached':'false'} if e is None else
                    {'IsAttached':'true','AttachedToObjID':'id-'+e[2],'AttachedToConnIndex':str(e[3])})
    with zipfile.ZipFile(path,'w') as z:z.writestr('saved.xml',ET.tostring(root))

def test_saved_topology_rejects_swapped_splitter_ports(tmp_path):
    path=tmp_path/'case.dwxmz';archive(path,CONNECTIONS)
    assert len(topology(path))==12
    changed={(a,1-i if a=='S-1' else i,b,j) for a,i,b,j in CONNECTIONS}
    archive(path,changed)
    with pytest.raises(RuntimeError,match='topology differs'):topology(path)

def test_altered_nested_process_artifact_rejected_before_replay(tmp_path):
    from coupling.run import sha
    f=tmp_path/'process/summary.json';f.parent.mkdir();f.write_text('{"power":1}\n')
    s={'status':'passed','artifact_sha256':{'process/summary.json':sha(f)}}
    (tmp_path/'summary.json').write_text(json.dumps(s));f.write_text('{"power":2}\n')
    with pytest.raises(RuntimeError,match='Artifact hash mismatch'):verify(tmp_path/'summary.json')

def test_actual_saved_dwsim_streams_match_real_julia_process_closure():
    # This compact record was exported from the actual solved/reloaded case;
    # fresh external solves are exercised by the documented study command.
    r=read(ROOT/'docs/process_filtration_reference.json');c=r['nominal_actual_DWSIM_case']
    source=r['source_for_Julia_replay'];p=r['protocol']
    _,closure=translate(runtime(),source,source['supports'][0],.5,20000.)
    evidence=validate(c['dwsim'],adapter_inputs(p,source,closure),closure,p,source)
    assert all(evidence['checks'].values())
    assert evidence['actual_specific_shaft_energy_kWh_m3_reference_volume']==pytest.approx(.02235249258998649)
    assert evidence['actual_recovery_mass_fraction']==pytest.approx(.5)

@pytest.mark.parametrize('stream,field,increment',[
    ('Permeate','mass_flow_kg_s',1e-6),('Retentate','enthalpy_J_kg',1.),
    ('Pump-discharge','viscosity_Pa_s',1e-4)])
def test_false_balance_flags_cannot_hide_changed_streams(stream,field,increment):
    r=read(ROOT/'docs/process_filtration_reference.json');c=r['nominal_actual_DWSIM_case']
    actual=copy.deepcopy(c['dwsim']);actual['streams'][stream][field]+=increment
    with pytest.raises((ValueError,RuntimeError)):
        validate(actual,actual['inputs'],c['julia_process'],r['protocol'],r['source_for_Julia_replay'])
