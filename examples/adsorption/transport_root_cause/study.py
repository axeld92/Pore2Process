"""Native matched-mesh runs with read-only peak and exact-population diagnostics."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from examples.adsorption.transport_accuracy.run import source_point, validate_record
from examples.adsorption.transport_accuracy.profile import analyze as profile_analysis

HERE = Path(__file__).parent
OUTPUT = ROOT/'results/adsorption/transport_root_cause'
DT = 1.9017455796408345e-5
END = .2147547125969688


def peak_analysis(path):
    path = Path(path)
    p = json.loads((path/'peak.json').read_text())
    f = np.genfromtxt(path/'peak.csv', delimiter=',', names=True)
    n = round(len(f)**(1/3))
    fluid = (f['material'] != 2).reshape((n,n,n))
    near = np.zeros_like(fluid)
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                if i or j or k:
                    near |= np.roll(~fluid, (i,j,k), axis=(0,1,2))
    near &= fluid
    c = f['concentration_kg_m3'].reshape(fluid.shape)
    amplitude = p['initial_amplitude_kg_m3']
    violation = ((c < .5-amplitude-.02) | (c > .5+amplitude+.02)) & fluid
    def location(index):
        row=f[index]
        return dict(ijk=[int(row[x]) for x in ('i','j','k')],
                    xyz_m=[float(row[x+'_m']) for x in ('x','y','z')],
                    concentration_kg_m3=float(row['concentration_kg_m3']),
                    touches_solid_link=bool(near.reshape(-1)[index]))
    ids=np.flatnonzero(fluid.reshape(-1))
    p.update(minimum=location(ids[np.argmin(c.reshape(-1)[ids])]),
             maximum=location(ids[np.argmax(c.reshape(-1)[ids])]),
             fraction_fluid_touching_solid_link=float(np.mean(near[fluid])),
             fraction_fluid_outside_engineering_bounds=float(np.mean(violation[fluid])),
             fraction_violators_touching_solid_link=None if not np.any(violation) else float(np.mean(near[violation])))
    return p


def analyze(path):
    path=Path(path)
    result=profile_analysis(path)
    result['peak']=peak_analysis(path)
    amplitude=result['peak']['initial_amplitude_kg_m3']
    result['checks']['lower_bound']=result['minimum_concentration_kg_m3']>=.5-amplitude-.02
    result['checks']['upper_bound']=result['maximum_concentration_kg_m3']<=.5+amplitude+.02
    result['passed']=all(result['checks'].values())
    result['constant_maximum_relative_departure']=None if amplitude else max(abs(result['minimum_concentration_kg_m3']/.5-1),abs(result['maximum_concentration_kg_m3']/.5-1))
    a=np.genfromtxt(path/'moments.csv', delimiter=',', names=True)
    a=np.atleast_1d(a)
    result['periodic_wave_final']={k:float(a[k][-1]) for k in a.dtype.names}
    return result


def execute(base, spec, binary):
    dest=base/'native'/spec['id'];dest.mkdir(parents=True)
    command=[str(binary),str(spec['n']),'1',str(spec['side_m']),'1e-9',
             str(.5+3e-9*DT/(spec['side_m']/spec['n'])**2), str(spec.get('end_s',END)),
             spec['field'],str(spec['speed_m_s']),'.5',str(spec['velocity_scale']),'1','zero',
             str(spec.get('audit_step',-1)),str(spec.get('amplitude',.25))]
    write_json(dest/'request.json',dict(spec=spec,command=command,binary_sha256=sha(binary),
               native_field_sha256=sha(Path(spec['field']))))
    started=perf_counter();print('Native start',spec['id'],flush=True)
    with (dest/'stdout.log').open('w') as log:
        p=subprocess.run(command,cwd=dest,stdout=log,stderr=subprocess.STDOUT,timeout=7200)
    write_json(dest/'execution.json',dict(returncode=p.returncode,elapsed_s=perf_counter()-started))
    if p.returncode:raise RuntimeError('Native failed: '+str(dest/'stdout.log'))
    print('Native finished',spec['id'],round(perf_counter()-started,2),'s',flush=True)


def check(path):
    path=Path(path).resolve();s=validate_record(path)
    for role in ('parent','hydraulic'):
        p=ROOT/s[role+'_summary']
        if sha(p)!=s[role+'_summary_sha256']:raise ValueError('Changed parent '+role)
        validate_record(p)
    for spec in s['specifications']:
        if sha(Path(spec['field']))!=spec['field_sha256']:raise ValueError('Changed velocity field')
        if analyze(path.parent/'native'/spec['id'])!=s['cases'][spec['id']]:raise ValueError('Analysis does not replay')
    print('Checked',len(s['source_sha256']),'sources,',len(s['artifacts_sha256']),'artifacts; 0 new native/DWSIM runs',flush=True)
    return s


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--grids',nargs='+',type=int,choices=[32,48,64,96],default=[32,48])
    p.add_argument('--point-json',type=Path)
    p.add_argument('--constant',action='store_true')
    p.add_argument('--audit-step',type=int,default=-1)
    p.add_argument('--end-s',type=float,default=END)
    p.add_argument('--check',type=Path)
    args=p.parse_args()
    if args.check:check(args.check);return
    parent_path=ROOT/json.loads((ROOT/'results/adsorption/transport_accuracy/diagnosis_latest.json').read_text())['summary']
    parent=validate_record(parent_path)
    hydraulic_path=ROOT/parent['hydraulic_summary'];validate_record(hydraulic_path)
    if sha(hydraulic_path)!=parent['hydraulic_summary_sha256']:raise ValueError('Hydraulic parent changed')
    supplied=None if args.point_json is None else json.loads(args.point_json.read_text())
    jobs=[]
    for n in args.grids:
        point=source_point(hydraulic_path.parent,'d800',n) if n in (32,48) else supplied
        if point is None or ('n' in point and point['n']!=n):raise ValueError('Matching new source point is required')
        if sha(Path(point['field']))!=point['field_sha256']:raise ValueError('Source point field hash mismatch')
        jobs.append(dict(point,id=f'n{n}'+('_constant' if args.constant else '_sine')+(f'_audit{args.audit_step}' if args.audit_step>=0 else ''),
                         n=n,amplitude=0. if args.constant else .25,audit_step=args.audit_step,end_s=args.end_s))
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    headers=[HERE/'flux_audit/population_audit.h',HERE/'flux_audit/uniform_advection_audit.h']
    files=[Path(__file__),HERE/'protocol.md',*headers,*list((HERE/'openlb').iterdir()),
           ROOT/'examples/adsorption/transport_accuracy/profile.py',ROOT/'examples/adsorption/transport_accuracy/run.py',
           ROOT/'examples/adsorption/column_transport/run.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    generated=base/'generated_source';shutil.copytree(HERE/'openlb',generated)
    for header in headers:shutil.copy2(header,generated/header.name)
    metadata=build_case(generated,ROOT/'.local/adsorption/transport_root_cause',base/'build')
    binary=ROOT/'.local/adsorption/transport_root_cause/capillary3d'
    record=dict(schema='pore2process.transport-root-cause.v1',run_id=rid,source_sha256=sources,
                specifications=jobs,build=metadata,parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
                hydraulic_summary=str(hydraulic_path.relative_to(ROOT)),hydraulic_summary_sha256=sha(hydraulic_path),
                new_native_runs=len(jobs),new_DWSIM_runs=0,new_adsorption_predictions=0,service_prediction_qualified=False)
    if args.point_json:
        shutil.copy2(args.point_json,base/'source_point.json')
        record['source_point_record']=str(args.point_json.resolve().relative_to(ROOT))
        record['source_point_record_sha256']=sha(args.point_json)
    write_json(base/'plan.json',record)
    print('Run',base,flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(lambda job:execute(base,job,binary),jobs))
    record['cases']={s['id']:analyze(base/'native'/s['id']) for s in jobs}
    record['vendor']=verify_vendor()
    record['artifacts_sha256']={str(x.relative_to(base)):sha(x) for x in base.rglob('*') if x.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Completed',base,flush=True)


if __name__=='__main__':main()
