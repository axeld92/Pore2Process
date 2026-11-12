"""Execute a bounded stock-equilibrium hypothesis test without editing parents."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from ...channel.run import PHYSICAL, analyze as analyze_channel
from ..analysis import summarize

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_root_cause/second_order_candidate'
PARENT=ROOT/'results/adsorption/transport_root_cause/runs/20261110T205424Z-ed34b563/native/n32_sine/request.json'


def analyze_porous(dest):
    raw=json.loads((dest/'raw.json').read_text())
    history=np.genfromtxt(dest/'inventory.csv',delimiter=',',names=True)
    request=json.loads((dest/'request.json').read_text())
    amplitude=request['amplitude_kg_m3']
    lo=float(np.min(history['min_concentration_kg_m3']));hi=float(np.max(history['max_concentration_kg_m3']))
    departure=max(.5-amplitude-lo,hi-(.5+amplitude),0.)
    initial=raw['initial_inventory_kg']
    drift=float(np.max(np.abs(history['mass_residual_kg']))/initial)
    audits={p.name:summarize(p,dx_m=raw['dx_m'],dt_s=raw['dt_s']) for p in dest.glob('audit_*.csv')}
    return dict(raw=raw,minimum_concentration_kg_m3=lo,maximum_concentration_kg_m3=hi,
                maximum_bound_excursion_kg_m3=departure,relative_inventory_drift=drift,
                checks=dict(finite=bool(np.isfinite(lo) and np.isfinite(hi)),mass=drift<=1e-8,
                            engineering_bounds=departure<=.02,
                            exact_population_balances=all(all(v['checks'].values()) for v in audits.values()),
                            **({'constant_preservation':departure<=1e-9} if amplitude==0 else {})),
                population_audits=audits,application_transport_qualified=False)


def run():
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    source_files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    source_files += [HERE.parent/'analysis.py',HERE.parent.parent/'channel/run.py',HERE.parent.parent/'channel/reference.py']
    sources=source_manifest(source_files)
    for name in sources:
        target=base/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    parent=json.loads(PARENT.read_text());shutil.copy2(PARENT,base/'parent_request.json')
    build={};binary={}
    for model in ('porous','channel'):
        builddir=ROOT/'.local/adsorption/transport_root_cause_secondorder'/model
        build[model]=build_case(HERE/(model+'_openlb'),builddir,base/'build'/model)
        binary[model]=builddir/'capillary3d'
    specifications=[('porous_sine','porous',.25,.75),('porous_constant','porous',0.,.75),
                    ('shear_sine','channel',.25,.75),('shear_constant','channel',0.,.75),
                    ('uniform_sine','channel',.25,0.)]
    cases={}
    for name,model,amplitude,shear in specifications:
        dest=base/'native'/name;dest.mkdir(parents=True)
        if model=='porous':
            command=[str(binary[model]),*parent['command'][1:]]
            command[-1]=str(amplitude)
            command[-2]='-1'
        else:
            p=PHYSICAL
            command=[str(binary[model]),'32',str(p['side_m']),str(p['diffusivity_m2_s']),str(p['dt_s']),
                     str(p['end_s']),str(p['speed_m_s']),str(shear),str(amplitude)]
        write_json(dest/'request.json',dict(command=command,model=model,amplitude_kg_m3=amplitude,
                   shear_fraction=shear,binary_sha256=sha(binary[model]),parent_request_sha256=sha(PARENT)))
        started=perf_counter();print('Second-order candidate starts',name,flush=True)
        with (dest/'stdout.log').open('w') as stream:
            result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=3600)
        write_json(dest/'execution.json',dict(returncode=result.returncode,wall_s=perf_counter()-started))
        if result.returncode:raise RuntimeError('Candidate native run failed: '+str(dest/'stdout.log'))
        cases[name]=analyze_porous(dest) if model=='porous' else analyze_channel(dest)
        write_json(dest/'analysis.json',cases[name]);print('Second-order candidate done',name,f'{perf_counter()-started:.1f}s',flush=True)
    for name,digest in sources.items():
        if sha(ROOT/name)!=digest:raise ValueError('Candidate source changed during run: '+name)
    record=dict(schema='pore2process.second-order-equilibrium-diagnostic.v1',run_id=rid,
                source_sha256=sources,build=build,vendor=verify_vendor(),cases=cases,
                actual_native_runs=len(cases),new_DWSIM_runs=0,application_transport_qualified=False,
                interpretation='Unqualified stock-composition hypothesis test; no numerical acceptance gate relaxed.')
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print(base/'summary.json',flush=True);return base/'summary.json'


def check():
    path=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];record=json.loads(path.read_text());base=path.parent
    for name,digest in record['source_sha256'].items():
        if sha(ROOT/name)!=digest or sha(base/'source_snapshot'/name)!=digest:raise ValueError('Candidate source changed: '+name)
    for name,digest in record['artifacts_sha256'].items():
        if sha(base/name)!=digest:raise ValueError('Candidate artifact changed: '+name)
    for name,expected in record['cases'].items():
        dest=base/'native'/name
        got=analyze_porous(dest) if name.startswith('porous') else analyze_channel(dest)
        if got!=expected:raise ValueError('Candidate analysis changed: '+name)
    result=dict(run_id=record['run_id'],sources_checked=len(record['source_sha256']),
                artifacts_checked=len(record['artifacts_sha256']),new_native_runs=0,new_DWSIM_runs=0,vendor=verify_vendor())
    write_json(OUTPUT/'verification.json',result);print(json.dumps(result,indent=2));return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    check() if args.check else run()


if __name__=='__main__':main()
