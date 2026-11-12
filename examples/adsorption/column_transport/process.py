"""Actual native -> Julia operating-point translation -> actual DWSIM pumps."""
from datetime import datetime,timezone
from functools import lru_cache
import argparse
import json
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4

import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.bridge import runtime,_check_thread,runtime_versions
from examples.adsorption.bridge import module as adsorption_module
from examples.adsorption.dwsim import build_adapter,solve_case
from .hydraulics import pump_inputs,OUTPUT,HERE
from .bridge import qualified_dispersion
from .check import check_hydraulics,check_hydraulic_diagnostics,integrity,load_latest


@lru_cache(maxsize=1)
def _loaded():
    adsorption_module()
    jl=runtime();jl.include(str(ROOT/'julia/adsorption/OperatingPointHydraulics.jl'))
    return jl.OperatingPointHydraulics


def translate(**inputs):
    _check_thread();jl=runtime();r=_loaded().translate(**inputs)
    return {str(k):float(getattr(r,str(k))) for k in jl.fieldnames(jl.typeof(r))}


def check():
    path=load_latest('process');s=json.loads(path.read_text());base=path.parent
    counts=integrity(s,base)
    native=check_hydraulics(ROOT/s['native_summary'])
    diagnostic=check_hydraulic_diagnostics()
    if diagnostic['source_summary_sha256']!=s['diagnostic_summary_sha256']:raise ValueError('Hydraulic diagnostic source changed.')
    if native['source_summary_sha256']!=s['native_summary_sha256']:raise ValueError('Native parent changed.')
    for label,c in s['cases'].items():
        replay=translate(**c['julia_inputs'])
        for key,value in replay.items():np.testing.assert_allclose(value,c['julia_output'][key],rtol=1e-12,atol=1e-20)
        pump=json.loads((base/label/'dwsim/summary.json').read_text())
        if pump!=c['dwsim']:raise ValueError('Actual DWSIM saved output changed.')
        np.testing.assert_allclose(pump['inputs']['Outlet_pressure_Pa']-pump['inputs']['Inlet_pressure_Pa'],replay['pressure_drop_Pa'],rtol=1e-10)
        np.testing.assert_allclose(pump['inputs']['Mass_flow_kg_s']/pump['inlet_density_kg_m3'],replay['flow_m3_s'],rtol=1e-12)
        np.testing.assert_allclose(pump['shaft_power_W'],replay['reference_shaft_power_W'],rtol=1e-8)
        assert all(c['checks'].values())
    assert not s['adsorption_qualified'] and s['service_time_s'] is None
    ref=dict(run_id=s['run_id'],source_summary=str(path.relative_to(ROOT)),source_summary_sha256=sha(path),
        native_summary=s['native_summary'],native_summary_sha256=s['native_summary_sha256'],
        cases={key:{k:c[k] for k in ('julia_inputs','julia_output','checks','shaft_power_W','specific_shaft_energy_Wh_m3')} for key,c in s['cases'].items()},
        sensitivities=s['sensitivities'],actual_new_DWSIM_cases=s['new_DWSIM_runs'],new_native_runs=s['new_native_runs'],
        new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,verified=counts,
        adsorption_qualified=False,service_time_s=None,service_volume_m3=None,
        native_to_Julia_to_DWSIM_executed=True)
    write_json(HERE/'reference/process_reference.json',ref)
    print('Actual Julia/DWSIM process checked',counts,'new simulator runs in check: 0',flush=True)
    return ref


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.check:check();return
    native_path=load_latest('hydraulics');native_ref=check_hydraulics(native_path)
    if not native_ref['hydraulic_checks_passed']:raise ValueError('Qualify native hydraulics first.')
    diagnostic=check_hydraulic_diagnostics()
    native=json.loads(native_path.read_text());p=native['inputs']
    closure=json.loads(load_latest('transport').read_text())['closure']
    try:qualified_dispersion(closure)
    except ValueError as e:reason=str(e)
    else:raise ValueError('This bounded process protocol must be updated before adding a new adsorption study.')
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'process/runs'/rid;base.mkdir(parents=True)
    files=[Path(__file__),HERE/'process_protocol.md',ROOT/'julia/adsorption/OperatingPointHydraulics.jl',
        ROOT/'julia/adsorption/test_operating_point_hydraulics.jl',ROOT/'julia/adsorption/AdsorptionColumn.jl',
        ROOT/'examples/adsorption/bridge.py',HERE/'hydraulics.py',HERE/'bridge.py',HERE/'check.py',
        HERE/'hydraulic_diagnostic_protocol.md',HERE/'hydraulic_diagnostics.py',
        *list((HERE/'openlb_domain').iterdir()),
        ROOT/'examples/adsorption/dwsim.py',*list((ROOT/'examples/adsorption/dwsim_runner').iterdir())]
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    build=build_adapter(base/'adapter_build');cold=perf_counter();_loaded();cold=perf_counter()-cold
    cases={};sensitivities=[]
    anchor=native['cases']['d800']['dwsim']
    for label,n in native['cases'].items():
        a=n['native']['48']['accepted']['analysis'];r=a['raw']
        parameters=dict(native_gradient_Pa_m=r['gradient_Pa_m'],native_velocity_m_s=a['superficial_velocity_m_s'],
            native_viscosity_Pa_s=r['viscosity_Pa_s'],native_density_kg_m3=r['density_kg_m3'],
            carrier_viscosity_Pa_s=p['carrier_viscosity_Pa_s'],carrier_density_kg_m3=p['carrier_density_kg_m3'],
            superficial_velocity_m_s=p['superficial_velocity_m_s'],length_m=p['length_m'],area_m2=p['area_m2'],efficiency_fraction=.75)
        point=translate(**parameters)
        inputs=pump_inputs(anchor,point['pressure_drop_Pa'],point['flow_m3_s'])
        inputs['Scope']='Actual native FCC operating-point resistance -> Julia load -> actual DWSIM bed-friction pump. No qualified new adsorption/service prediction; zero elevation and ancillary duty.'
        pump=solve_case(base/label/'dwsim',inputs)
        checks=dict(work=abs(pump['shaft_power_W']/point['reference_shaft_power_W']-1)<1e-8,
            mass=abs(pump['mass_balance_error_kg_s'])<=1e-10*pump['inlet_mass_flow_kg_s'],
            enthalpy=abs(pump['energy_balance_error_W'])<=.001*pump['shaft_power_W'],
            reload=pump['save_reload_solve_checked'],
            native_Reynolds_match=point['Reynolds_relative_mismatch']<=1e-4)
        cases[label]=dict(julia_inputs=parameters,julia_output=point,dwsim=pump,checks=checks,
            shaft_power_W=pump['shaft_power_W'],specific_shaft_energy_Wh_m3=pump['shaft_power_W']/point['flow_m3_s']/3600)
        for efficiency in (.6,.75,.85):
            for extra in (0.,1000.,10000.):
                sensitivities.append(dict(particle=label,efficiency=efficiency,extra_pressure_Pa=extra,
                    shaft_energy_Wh_m3=(point['pressure_drop_Pa']+extra)/efficiency/3600,
                    role='Algebraic chosen scenario using verified pump-work law; not another DWSIM solve or uncertainty interval.'))
    s=dict(schema='pore2process.native-julia-dwsim-point.v1',run_id=rid,source_sha256=sources,
        native_summary=str(native_path.relative_to(ROOT)),native_summary_sha256=sha(native_path),
        diagnostic_summary=diagnostic['source_summary'],diagnostic_summary_sha256=diagnostic['source_summary_sha256'],
        cases=cases,checks_passed=all(all(c['checks'].values()) for c in cases.values()),sensitivities=sensitivities,
        new_native_runs=0,new_DWSIM_runs=2,new_adsorption_study_runs=0,
        adsorption_qualified=False,adsorption_blocker=reason,service_time_s=None,service_volume_m3=None,
        build=build,versions=runtime_versions(),runtime_and_module_initialization_s=cold)
    s['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',s)
    write_json(OUTPUT/'process/latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Actual native -> Julia -> DWSIM process saved',base,s['checks_passed'],flush=True)


if __name__=='__main__':main()
