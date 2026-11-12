"""Saved-evidence checks. Replays Julia analysis, not native/DWSIM simulations."""
import json
from pathlib import Path

import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import verified_source_path
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.bed import analyze as flow_analysis
from .bridge import qualified_dispersion
from .run import analyze, OUTPUT


def load_latest(stage):
    pointer=OUTPUT/stage/'latest.json'
    return ROOT/json.loads(pointer.read_text())['summary']


def integrity(summary,base):
    for name,digest in summary['source_sha256'].items():
        source=verified_source_path(name,digest)
        if sha(source)!=digest:raise ValueError('Implementation changed: '+name)
    for name,digest in summary['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Saved artifact changed: '+name)
    return dict(source_files=len(summary['source_sha256']),artifacts=len(summary['artifacts_sha256']))


def check_transport(path=None):
    path=Path(path) if path else load_latest('transport')
    s=json.loads(path.read_text());base=path.parent
    checked=integrity(s,base)
    for key,saved in s['cases'].items():
        actual=analyze(base/'native'/key)
        if actual!=saved:raise ValueError('Native analysis differs: '+key)
    failures=[k for k,c in s['cases'].items() if not c['qualification_passed']]
    try:
        qualified_dispersion(s['closure'])
    except ValueError:
        rejected=True
    else:
        rejected=False
    if failures and (s['transport_qualified'] or not rejected):
        raise ValueError('Failed native transport was promoted to adsorption.')
    if s['closure']['axial_dispersion_m2_s'] is not None and failures:
        raise ValueError('A coefficient was extracted from failed native transport.')
    compact=dict(schema=s['schema'],run_id=s['run_id'],source_summary=str(path.relative_to(ROOT)),
        source_summary_sha256=sha(path),transport_qualified=s['transport_qualified'],
        application=s['application'],closure=s['closure'],native_runs=s['native_runs'],
        cases={k:{name:c[name] for name in ('checks','minimum_concentration_kg_m3',
            'maximum_concentration_kg_m3','uniform_maximum_departure','normalized_mass_residual','qualification_passed')}
            | ({'comparison':c['comparison']} if 'comparison' in c else {}) for k,c in s['cases'].items()},
        failed_native_cases=failures,unqualified_handoff_rejected=rejected,
        new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,verified=checked)
    write_json(Path(__file__).parent/'reference/transport_reference.json',compact)
    return compact


def check_hydraulics(path=None):
    path=Path(path) if path else load_latest('hydraulics')
    s=json.loads(path.read_text());base=path.parent
    checked=integrity(s,base);p=s['inputs']
    mass=(1-p['nominal_porosity'])*p['particle_density_kg_m3']*p['area_m2']*p['length_m']
    np.testing.assert_allclose(mass,p['dry_mass_kg'],rtol=1e-13)
    rows={}
    for key,c in s['cases'].items():
        if not all(c['checks'].values()):raise ValueError('Hydraulic qualification failed: '+key)
        for group in c['native'].values():
            for iteration in group['iterations']:
                if flow_analysis(base/iteration['directory'])!=iteration['analysis']:
                    raise ValueError('Native hydraulic field/history changed.')
        a=c['native']['48']['accepted']['analysis']
        g=(a['raw']['gradient_Pa_m']/a['superficial_velocity_m_s'])*p['native_target_velocity_m_s']
        nu=p['carrier_viscosity_Pa_s']/p['carrier_density_kg_m3']
        pressure=g*(p['carrier_density_kg_m3']/1000)*(nu/1e-6)**2*p['length_m']
        np.testing.assert_allclose(pressure,c['pressure_drop_Pa'],rtol=1e-13)
        pump=json.loads((base/key/'dwsim/summary.json').read_text())
        if pump!=c['dwsim']:raise ValueError('DWSIM saved summary differs.')
        actual_flow=pump['inlet_mass_flow_kg_s']/pump['inlet_density_kg_m3']
        expected=actual_flow*pressure/.75
        np.testing.assert_allclose(actual_flow,p['flow_m3_s'],rtol=1e-12)
        np.testing.assert_allclose(pump['outlet_pressure_Pa']-pump['inlet_pressure_Pa'],pressure,rtol=1e-10)
        np.testing.assert_allclose(pump['shaft_power_W'],expected,rtol=1e-8)
        assert pump['save_reload_solve_checked']
        if any(c[k] is not None for k in ('breakthrough_time_s','service_volume_m3','cumulative_service_energy_kWh')):
            raise ValueError('Unqualified service prediction must remain unavailable.')
        rows[key]={k:c[k] for k in ('diameter_m','particle_Reynolds','pressure_drop_Pa','pressure_gradient_Pa_m',
            'shaft_power_W','specific_pump_energy_kWh_m3','refinement','checks','breakthrough_time_s',
            'service_volume_m3','cumulative_service_energy_kWh')}
    ratio=rows['d450']['shaft_power_W']/rows['d800']['shaft_power_W']
    compact=dict(schema=s['schema'],run_id=s['run_id'],source_summary=str(path.relative_to(ROOT)),
        source_summary_sha256=sha(path),inputs=p,cases=rows,power_ratio_small_to_large=ratio,
        hydraulic_checks_passed=s['hydraulic_checks_passed'],service_prediction_qualified=False,
        actual_native_flow_runs=s['new_native_flow_runs']+s.get('reused_native_flow_runs',0),
        actual_DWSIM_cases=s['new_DWSIM_runs']+s.get('reused_DWSIM_cases',0),
        new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,verified=checked)
    write_json(Path(__file__).parent/'reference/hydraulic_reference.json',compact)
    return compact


def check_sensitivity():
    path=load_latest('sensitivity');s=json.loads(path.read_text());counts=integrity(s,path.parent)
    for key,value in s['cases'].items():
        if analyze(path.parent/'native'/key)!=value:raise ValueError('Preservation sensitivity changed.')
    compact=dict(run_id=s['run_id'],source_summary=str(path.relative_to(ROOT)),source_summary_sha256=sha(path),
        original_departure=s['original_departure'],preserved=s['preserved'],native_runs=s['native_runs'],
        cases={k:{q:c[q] for q in ('uniform_maximum_departure','minimum_concentration_kg_m3',
            'maximum_concentration_kg_m3','normalized_mass_residual','checks')} for k,c in s['cases'].items()},
        verified=counts,new_native_runs_in_check=0,new_DWSIM_runs_in_check=0)
    write_json(Path(__file__).parent/'reference/sensitivity_reference.json',compact)
    return compact


def check_hydraulic_diagnostics():
    path=load_latest('hydraulic_diagnostics');s=json.loads(path.read_text());counts=integrity(s,path.parent)
    reference=flow_analysis(ROOT/s['source_reference'])
    if sha(ROOT/s['source_reference']/'field.csv')!=s['reference_field_sha256']:
        raise ValueError('Hydraulic initial/domain reference changed.')
    for name,c in s['cases'].items():
        raw=json.loads((path.parent/name/'raw.json').read_text())
        a=np.genfromtxt(path.parent/name/'field.csv',delimiter=',',names=True)
        n=raw['resolution'];cells=raw['domain_cells_x']
        u=float(np.sum(a['ux_m_s'])/(cells*n**3))
        np.testing.assert_allclose(u,c['superficial_velocity_m_s'],rtol=1e-13)
        delta=abs(u/reference['superficial_velocity_m_s']-1)
        np.testing.assert_allclose(delta,c['mean_velocity_relative_difference'],rtol=1e-10,atol=1e-15)
        if cells==2:
            values=np.zeros((2*n,n,n,3))
            i,j,k=(a[q].astype(int) for q in ('i','j','k'))
            for axis,q in enumerate(('ux_m_s','uy_m_s','uz_m_s')):values[i,j,k,axis]=a[q]
            independent=float(np.sqrt(np.sum((values[:n]-values[n:])**2)/np.sum(values[:n]**2)))
            np.testing.assert_allclose(independent,raw['periodicity_relative_error'],rtol=1e-8,atol=1e-14)
        if not all(c['checks'].values()):
            raise ValueError('Hydraulic initial/domain qualification failed: '+name)
    ref=dict(run_id=s['run_id'],source_summary=str(path.relative_to(ROOT)),source_summary_sha256=sha(path),
        checks_passed=s['checks_passed'],native_runs=s['native_runs']+s.get('reused_native_runs',0),cases=s['cases'],verified=counts,
        new_native_runs_in_check=0,new_DWSIM_runs_in_check=0)
    write_json(Path(__file__).parent/'reference/hydraulic_diagnostic_reference.json',ref)
    return ref


def main():
    t=check_transport()
    h=check_hydraulics()
    sensitivity=check_sensitivity()
    diagnostics=check_hydraulic_diagnostics()
    process=None
    if (OUTPUT/'process/latest.json').exists():
        from .process import check as process_check
        process=process_check()
    print(json.dumps(dict(saved_transport=t['verified'],saved_hydraulics=h['verified'],
        saved_sensitivity=sensitivity['verified'],saved_process=None if process is None else process['verified'],
        saved_hydraulic_diagnostics=diagnostics['verified'],
        transport_qualified=t['transport_qualified'],hydraulics_qualified=h['hydraulic_checks_passed'],
        service_prediction_qualified=False,new_native_runs=0,new_DWSIM_runs=0,
        vendor=verify_vendor()),indent=2))


if __name__=='__main__':main()
