"""Execute and independently verify stock-native periodic shear controls."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from .reference import converged_reference

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_root_cause/channel'
PHYSICAL=dict(side_m=.0013891725543141453,speed_m_s=.012937295182166382,
              diffusivity_m2_s=1e-9,dt_s=1.90174557964e-5,end_s=.214764128309)


def specifications():
    return [dict(id='uniform32',n=32,shear_fraction=0.,amplitude=.25),
            dict(id='constant_shear32',n=32,shear_fraction=.75,amplitude=0.),
            *[dict(id=f'shear{n}',n=n,shear_fraction=.75,amplitude=.25) for n in (32,48,64)]]


def execute(base,spec,binary):
    dest=base/'native'/spec['id'];dest.mkdir(parents=True)
    p=PHYSICAL
    command=[str(binary),str(spec['n']),str(p['side_m']),str(p['diffusivity_m2_s']),str(p['dt_s']),
             str(p['end_s']),str(p['speed_m_s']),str(spec['shear_fraction']),str(spec['amplitude'])]
    write_json(dest/'request.json',dict(spec=spec,physical=p,command=command,binary_sha256=sha(binary)))
    print('Native shear start',spec['id'],flush=True);started=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=3600)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-started))
    if result.returncode:raise RuntimeError('Native shear run failed: '+str(dest/'stdout.log'))
    print('Native shear finished',spec['id'],f'{perf_counter()-started:.1f}s',flush=True)


def analyze(dest):
    dest=Path(dest);raw=json.loads((dest/'raw.json').read_text())
    a=np.genfromtxt(dest/'inventory.csv',delimiter=',',names=True)
    f=np.genfromtxt(dest/'final.csv',delimiter=',',names=True)
    low=float(np.min(a['min_concentration_kg_m3']));high=float(np.max(a['max_concentration_kg_m3']))
    amplitude=raw['initial_amplitude_kg_m3'];mass=float(np.max(abs(a['inventory_kg']/a['inventory_kg'][0]-1)))
    departure=max(.5-amplitude-low,high-(.5+amplitude),0.)
    checks=dict(finite=bool(all(np.all(np.isfinite(a[k])) for k in a.dtype.names)),
                mass=mass<=1e-8,z_invariance=float(np.max(f['z_variation_kg_m3']))<=1e-12,
                strict_bounds=departure<=1e-8,engineering_bounds=departure<=.02)
    result=dict(raw=raw,checks=checks,minimum_concentration_kg_m3=low,maximum_concentration_kg_m3=high,
                maximum_bound_excursion_kg_m3=departure,relative_inventory_drift=mass,
                application_transport_qualified=False,new_service_prediction=False)
    if amplitude:
        physical=dict(side_m=raw['side_m'],speed_m_s=raw['speed_m_s'],shear_fraction=raw['shear_fraction'],
                      diffusivity_m2_s=raw['diffusivity_m2_s'],duration_s=raw['duration_s'])
        exact,reference=converged_reference(f['x_m'],f['y_m'],amplitude_kg_m3=amplitude,**physical)
        difference=f['concentration_kg_m3']-exact
        result.update(reference_convergence=reference,RMS_error_over_amplitude=float(np.sqrt(np.mean(difference*difference)))/amplitude,
                      maximum_error_over_amplitude=float(np.max(abs(difference)))/amplitude)
        checks['profile']=result['RMS_error_over_amplitude']<=.01
        if not raw['shear_fraction']:
            k=2*np.pi/raw['side_m'];uniform=.5+amplitude*np.exp(-raw['diffusivity_m2_s']*k*k*raw['duration_s'])*np.sin(k*(f['x_m']-raw['speed_m_s']*raw['duration_s']))
            result['independent_uniform_reference_max_error_kg_m3']=float(np.max(abs(exact-uniform)))
            checks['uniform_oracle']=result['independent_uniform_reference_max_error_kg_m3']<=1e-10
    else:
        result['maximum_constant_error_kg_m3']=max(abs(low-.5),abs(high-.5))
        checks['constant_preservation']=result['maximum_constant_error_kg_m3']<=1e-9
    result['control_checks_passed']=all(checks.values())
    return result


def check(path=None):
    p=Path(path) if path else ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary']
    s=json.loads(p.read_text());base=p.parent
    for name,digest in s['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Current source changed: '+name)
        if sha(base/'source_snapshot'/name)!=digest:raise ValueError('Snapshot changed: '+name)
    for name,digest in s['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Artifact changed: '+name)
    for key,value in s['cases'].items():
        if analyze(base/'native'/key)!=value:raise ValueError('Analysis changed: '+key)
    vendor=verify_vendor()
    result=dict(run_id=s['run_id'],source_summary=str(p.relative_to(ROOT)),source_summary_sha256=sha(p),
                sources_checked=len(s['source_sha256']),artifacts_checked=len(s['artifacts_sha256']),
                actual_native_runs=len(s['cases']),new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,
                vendor=vendor,cases=s['cases'])
    write_json(OUTPUT/'verification.json',result);print(json.dumps(result,indent=2));return result


def run(workers=1):
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    sources=source_manifest([HERE/'run.py',HERE/'reference.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir())])
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    builddir=ROOT/'.local/adsorption/transport_root_cause_channel'
    metadata=build_case(HERE/'openlb',builddir,base/'build')
    jobs=specifications()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(lambda spec:execute(base,spec,builddir/'capillary3d'),jobs))
    record=dict(schema='pore2process.periodic-shear-diagnosis.v1',run_id=rid,source_sha256=sources,
                specifications=jobs,physical_inputs=PHYSICAL,build=metadata,vendor=verify_vendor(),
                cases={spec['id']:analyze(base/'native'/spec['id']) for spec in jobs},
                new_native_runs=len(jobs),new_DWSIM_runs=0,application_transport_qualified=False)
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Shear study completed',base,flush=True);return base/'summary.json'


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--workers',type=int,default=1,choices=(1,2))
    args=p.parse_args();check() if args.check else run(args.workers)

if __name__=='__main__':main()
