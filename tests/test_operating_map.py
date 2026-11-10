"""Scientific transfer bounds and real native/process contracts."""
import copy
from pathlib import Path
import numpy as np
import pytest
from pore2process_demo.bridge import runtime
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from examples.filtration.operating_map.evidence import extract
from examples.filtration.operating_map.model import capacity_decisions, translate
from examples.filtration.calibration.check import close_tree


def test_source_conditions_reading_bounds_and_no_pooling():
    e=extract();points=e['points']
    assert [x['source_Re'] for x in points]==[1388,2776,3701,4626]
    assert points[0]['critical_flux_LMH']==21.4
    assert points[2]['critical_flux_LMH']==72.5
    assert all(x['reading_interval_LMH'][0]<=x['critical_flux_LMH']<=x['reading_interval_LMH'][1] for x in points)
    assert all(x['experimental_point_uncertainty'] is None for x in points)
    assert all(x is None for x in e['unavailable'].values())
    assert len({s['doi'] for s in e['sources']})==4


def test_capacity_tests_peak_flux_and_keeps_extrapolation_unavailable():
    p=read(ROOT/'examples/filtration/operating_map/protocol.json');e=extract();jl=runtime()
    a={'axial':{'inlet_Re':1388.,'inlet_velocity_m_s':.5,'local_flux_m_s':[24/3.6e6,12/3.6e6]}}
    rows=capacity_decisions(jl,a,e,p)
    nominal=next(v for v in rows if (v['mapping'],v['capacity_fraction'],v['transfer_factor'])==('equal_Re',1.,1.))
    assert not nominal['nominal_admissible']
    assert nominal['maximum_local_flux_LMH']==pytest.approx(24)
    speed=[v for v in rows if v['mapping']=='equal_velocity']
    assert all(v['status']=='unavailable' and v['nominal_admissible'] is None for v in speed)
    assert nominal['sustainable_operation_validated'] is False


def test_independent_source_interpolation_and_shape_bounds():
    j=runtime();v=j.critical_flux_capacity([1000.,2000.],[20/3.6e6,50/3.6e6],
        [18/3.6e6,45/3.6e6],[22/3.6e6,55/3.6e6],1500.,30/3.6e6,.95,1.)
    assert float(v.nominal_capacity_m_s)==pytest.approx(.95*35/3.6e6)
    assert bool(v.nominal_admissible)
    assert str(v.monotonic_classification)=='indeterminate'


def reference():return read(ROOT/'docs/operating_map_reference.json')


def test_real_native_measurement_drives_saved_dwsim_replay():
    r=reference();c=r['actual_case'];q=c['request'];jl=runtime()
    t=translate(jl,r['native_for_replay'],r['TRT_support_for_replay'],r['protocol'],q)
    close_tree(t,c['translation'],'Real Julia replay')
    from examples.filtration.crossflow_recycle.model import validate
    close_tree(validate(c['dwsim'],c['dwsim']['inputs'],t['process'],r['hydraulics'],r['TRT_support_for_replay']),c['validation'])
    assert t['axial']['geometry_coefficient_m4']==r['native_for_replay']['geometry_coefficient_m4']
    assert t['axial']['inlet_Re']>1000
    assert t['sustainable_flux_LMH'] is None and t['fouling_time_s'] is None
    assert t['axial']['outlet_TMP_Pa']>0
    assert len(t['axial']['x_m'])==101
    assert c['dwsim']['inputs']['Force_liquid'] is True
    assert c['dwsim']['recycle']['maximum_iterations']==20000
    assert [v['assumed_local_rejection'] for v in t['quality_sensitivity']]==[.5,.8]
    assert t['quality_sensitivity'][0]['balance']['permeate_COD_kg_m3']>t['feed']['permeate_COD_kg_m3']>t['quality_sensitivity'][1]['balance']['permeate_COD_kg_m3']
    assert np.all(np.diff(t['axial']['pressure_Pa'])<0)


@pytest.mark.parametrize('change',[{'inlet_Re':1.},{'inlet_Re':3000.},{'flux_LMH':float('nan')},
    {'flux_LMH':100.},{'viscosity_multiplier':0.},{'fixed_multiplier':0.},{'support':'unknown'}])
def test_unsupported_physics_rejected_before_dwsim(change):
    r=reference();q={**r['actual_case']['request'],**change}
    with pytest.raises(ValueError):translate(runtime(),r['native_for_replay'],r['TRT_support_for_replay'],r['protocol'],q)


def test_changed_native_artifact_is_rejected(tmp_path):
    from pore2process_demo.io import sha
    from examples.filtration.operating_map.check import artifacts
    f=tmp_path/'native.csv';f.write_text('x,u\n1,2\n')
    s={'artifacts_sha256':{'native.csv':sha(f)}};f.write_text('x,u\n1,3\n')
    with pytest.raises(ValueError,match='Changed artifact'):artifacts(tmp_path/'summary.json',s)


def test_retained_measurements_replay_without_licensed_originals(monkeypatch,tmp_path):
    import examples.filtration.calibration.audit_crossflow_sources as audit_module
    before=extract();monkeypatch.setattr(audit_module,'ROOT',tmp_path)
    after=extract()
    assert before['points']==after['points']
    assert all(v['status'].startswith('unavailable_locally') for v in after['source_audit'])
    with pytest.raises(FileNotFoundError):extract(require_sources=True)
