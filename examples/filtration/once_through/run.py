"""Verified stored OpenLB designs → Julia loads → explicit DWSIM filtration."""
from pore2process_demo.paths import project_path
from pore2process_demo.provenance import implementation_files, source_manifest
from pore2process_demo.io import read_json as read
import argparse
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
from time import perf_counter
from uuid import uuid4

from pore2process_demo.io import sha
from examples.filtration.graded_support.check import verify as verify_source
from pore2process_demo.io import write_json
from pore2process_demo.bridge import runtime,runtime_versions
from examples.filtration.once_through import dwsim
from examples.filtration.once_through.model import adapter_inputs, translate, validate, topology

from pore2process_demo.paths import ROOT
def check_protocol(p,source):
    if p['protocol_id']!='once-through-filtration-v1' or source['protocol']['protocol_id']!=p['source_protocol_id']:
        raise RuntimeError('Changed study/source protocol; review before running')
    q=source['protocol']['process']
    for a,b in (('temperature_K','temperature_K'),('external_absolute_pressure_Pa','inlet_absolute_pressure_Pa'),
                ('pump_efficiency_fraction','pump_efficiency_fraction')):
        if p[a]!=q[b]: raise RuntimeError('Do not reuse native calibration with changed anchor conditions')
    if q['filtration_flux_LMH']!=100 or q['membrane_area_m2']!=1:
        raise RuntimeError('This scenario requires its declared fixed 100 LMH and 1 m2')

def jobs(p,source):
    for channel in p['channel_pressure_drops_Pa']:
        for recovery in p['permeate_recoveries']:
            for support in source['supports']:
                yield f"channel_{channel:g}Pa/recovery_{100*recovery:g}/{support['name']}",support,recovery,channel

def comparison_rows(cases):
    rows=[]
    for c in cases:
        j=c['julia_process'];d=c['validation']
        rows.append({'case_id':c['case_id'],'support':c['support'],'recovery_fraction':c['recovery_fraction'],
            'channel_pressure_drop_Pa':c['channel_pressure_drop_Pa'],
            'transmembrane_pressure_Pa':j['transmembrane_pressure_Pa'],
            'permeate_reference_flow_m3_s':j['permeate_reference_flow_m3_s'],
            'feed_reference_flow_m3_s':j['feed_reference_flow_m3_s'],
            'retentate_reference_flow_m3_s':j['retentate_reference_flow_m3_s'],
            'actual_DWSIM_shaft_power_W':c['dwsim']['shaft_power_W'],
            'actual_specific_shaft_energy_kWh_m3_reference_volume':d['actual_specific_shaft_energy_kWh_m3_reference_volume'],
            'energy_saving_vs_straight_fraction':c['energy_saving_vs_straight_fraction'],
            'max_mass_relative_residual':d['max_mass_relative_residual'],
            'max_energy_relative_residual':d['max_energy_relative_residual']})
    return rows

def execute(output,source_path):
    p=read(ROOT/'examples/filtration/once_through/protocol.json');source_path=Path(source_path).resolve();source=read(source_path)
    check_protocol(p,source)
    # Start and time the one shared runtime before the upstream checker uses it.
    t=perf_counter();jl=runtime();cold=perf_counter()-t
    print('Checking immutable native measurements, literature boundary and earlier process evidence...',flush=True)
    before=sha(source_path);verification=verify_source(source_path,write_record=False)
    if sha(source_path)!=before: raise RuntimeError('Source summary changed during verification')
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=Path(output)/'runs'/run_id;base.mkdir(parents=True)
    write_json(base/'status.json',{'status':'running','run_id':run_id})
    try:
        shutil.copy2(source_path,base/'source_summary.json')
        write_json(base/'source_verification.json',verification)
        adapter=dwsim.build(base/'adapter_build')
        implementation=sorted([f for f in (ROOT/'examples/filtration/once_through').rglob('*') if 'reference' not in f.parts and f.suffix in ('.py','.cs','.csproj','.json')])
        implementation += [project_path(ROOT/f) for f in ('examples/filtration/graded_support/check.py','examples/filtration/graded_support/model.py','julia/src/ProcessFiltration.jl',
            'julia/src/GradedSupport.jl','julia/src/SeriesCalibration.jl','julia/src/Pore2ProcessDemo.jl',
            'src/pore2process_demo/runtime/dwsim.py','src/pore2process_demo/bridge.py','uv.lock','julia/Project.toml','julia/Manifest.toml')]
        implementation = implementation_files(implementation)
        snapshot={'protocol':p,'source_summary_path':str(source_path.relative_to(ROOT)),
            'source_summary_sha256':before,'source_run_id':source['run_id'],
            'implementation_sha256':source_manifest(implementation),
            'runtime_sha256':{str(f.relative_to(ROOT)):sha(f) for f in (dwsim.DISTRIBUTION/name for name in
                ('DWSIM.Automation.FluentAPI.dll','DWSIM.UnitOperations.dll','DWSIM.Thermodynamics.dll','DWSIM.Interfaces.dll'))},
            'machine':{'system':platform.system(),'release':platform.release(),'architecture':platform.machine()}}
        write_json(base/'inputs.json',snapshot)
        cases=[]
        for name,support,recovery,channel in jobs(p,source):
            load,closure=translate(jl,source,support,recovery,channel)
            inputs=adapter_inputs(p,source,closure)
            print(f'Solving/reloading actual DWSIM filtration: {name}',flush=True)
            directory=base/'process'/name;actual=dwsim.solve(directory,inputs)
            evidence=validate(actual,inputs,closure,p,source)
            connections=topology(directory/'examples.filtration.layered_support.dwxmz')
            baseline=next((c for c in cases if c['channel_pressure_drop_Pa']==channel and c['recovery_fraction']==recovery),None)
            saving=0. if baseline is None else 1-actual['shaft_power_W']/baseline['dwsim']['shaft_power_W']
            cases.append({'case_id':name,'support':support['name'],'recovery_fraction':recovery,
                'channel_pressure_drop_Pa':channel,'julia_membrane':load,'julia_process':closure,
                'dwsim':actual,'validation':evidence,'saved_connections':connections,
                'energy_saving_vs_straight_fraction':saving})
        rows=comparison_rows(cases)
        with (base/'comparison.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        from examples.filtration.once_through.plot import plot
        plot(cases,base/'energy_comparison.png',p['nominal_recovery'])
        nominal=base/'process/channel_20000Pa/recovery_50/straight/flowsheet.png'
        shutil.copy2(nominal,base/'flowsheet.png')
        summary={'schema_version':1,'status':'passed','run_id':run_id,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
            **snapshot,'source_verification':verification,'adapter':adapter,
            'versions':runtime_versions(),'Julia_initialization_s':cold,'cases':cases,
            'scope':'Explicit once-through pure-water hydraulic filtration; declared recovery, no recycle or pressure recovery; not total plant energy.',
            'artifact_sha256':{str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file() and f.name!='status.json'}}
        write_json(base/'summary.json',summary)
        from examples.filtration.once_through.check import verify
        verified=verify(base/'summary.json',verify_upstream=False)
        write_json(base/'verification.json',verified)
        write_json(base/'status.json',{'status':'passed','run_id':run_id})
        write_json(Path(output)/'latest.json',{'summary':str((base/'summary.json').resolve())})
        print(f'Passed: {base / "summary.json"}',flush=True)
        return summary
    except BaseException as e:
        write_json(base/'status.json',{'status':'failed','run_id':run_id,'error':str(e)})
        raise

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--output',type=Path,default=ROOT/'results/process_filtration')
    a.add_argument('--source',type=Path)
    args=a.parse_args()
    source=args.source or Path(read(ROOT/'results/graded_support/latest.json')['summary'])
    execute(args.output,source)

if __name__=='__main__': main()
