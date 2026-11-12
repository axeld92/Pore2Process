"""Fresh native time/mesh probes of the best failed preservation candidate."""
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
from uuid import uuid4
import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from examples.adsorption.bed import analyze as flow_analysis
from .run import execute,analyze,OUTPUT,HERE
from .check import load_latest


def main():
    prior_path=load_latest('transport');prior=json.loads(prior_path.read_text())
    request=prior_path.parent/'native/application_constant_regularized27/request.json'
    original=json.loads(request.read_text())
    binary=ROOT/'.local/adsorption/column_transport_model1/capillary3d'
    if sha(binary)!=original['binary_sha256']:raise ValueError('Native binary changed.')
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'sensitivity/runs'/rid;base.mkdir(parents=True)
    sources=source_manifest([Path(__file__),HERE/'sensitivity_protocol.md',HERE/'run.py',
        *list((HERE/'openlb').iterdir())])
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    jobs=[]
    for name,factor in [('half_dt',.5),('tenth_dt',.1)]:
        spec=dict(original['spec']);spec['id']=name
        spec['tau_diffusive']=.5+(spec['tau_diffusive']-.5)*factor
        jobs.append(spec)
    root=ROOT/prior['application']['hydraulic_run']
    candidates=sorted((root/'native_flow/d800_grid48').glob('iteration*/raw.json'),
        key=lambda p:int(p.parent.name.removeprefix('iteration')))
    field=candidates[-1].parent/'field.csv'
    flow=flow_analysis(field.parent)
    carrier=json.loads((root/'inputs.json').read_text())
    if abs(flow['superficial_velocity_m_s']/carrier['native_target_velocity_m_s']-1)>1e-4:
        raise ValueError('Need completed fixed-flow N48 field.')
    spec=dict(original['spec']);spec.update(id='grid48',n=48,field=str(field))
    h=spec['side_m']/48
    a=np.genfromtxt(field,delimiter=',',names=True)
    maxu=float(np.max(np.sqrt(sum(a[k]**2 for k in ('ux_m_s','uy_m_s','uz_m_s')))))*spec['velocity_scale']
    dt=.02*h/maxu
    spec['tau_diffusive']=.5+spec['D_m2_s']*dt*3/h**2
    jobs.append(spec)
    cases={}
    for spec in jobs:
        execute(base,spec,binary);cases[spec['id']]=analyze(base/'native'/spec['id'])
    result=dict(schema='pore2process.column-transport-sensitivity.v1',run_id=rid,
        source_sha256=sources,prior_summary_sha256=sha(prior_path),cases=cases,
        original_departure=prior['cases']['application_constant_regularized27']['uniform_maximum_departure'],
        preserved=all(c['checks']['uniform_preservation'] for c in cases.values()),native_runs=len(jobs),
        new_adsorption_predictions=0,new_DWSIM_runs=0)
    result['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',result)
    write_json(OUTPUT/'sensitivity/latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Sensitivity saved',base,{k:c['uniform_maximum_departure'] for k,c in cases.items()},flush=True)


if __name__=='__main__':main()
