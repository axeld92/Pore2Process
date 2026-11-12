"""Actual cold/domain native tests of the higher-Re operating point."""
from datetime import datetime,timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from uuid import uuid4
import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from examples.adsorption.bed import analyze
from .hydraulics import OUTPUT,HERE


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--hydraulic-run',type=Path)
    parser.add_argument('--reuse-native',type=Path)
    args=parser.parse_args()
    if args.hydraulic_run:
        hydraulic=args.hydraulic_run.resolve()
    else:
        pointer=json.loads((OUTPUT/'hydraulics/latest.json').read_text())
        hydraulic=(ROOT/pointer['summary']).parent
    candidates=sorted((hydraulic/'native_flow/d800_grid32').glob('iteration*/raw.json'),
        key=lambda p:int(p.parent.name.removeprefix('iteration')))
    if not candidates:raise ValueError('Run the fixed-flow d800 N32 native group first.')
    reference=candidates[-1].parent
    old=analyze(reference);r=old['raw']
    target=json.loads((hydraulic/'inputs.json').read_text())['native_target_velocity_m_s']
    if abs(old['superficial_velocity_m_s']/target-1)>1e-4:
        raise ValueError('Reference native group has not reached the imposed flow.')
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'hydraulic_diagnostics/runs'/rid;base.mkdir(parents=True)
    sources=source_manifest([Path(__file__),HERE/'hydraulic_diagnostic_protocol.md',*list((HERE/'openlb_domain').iterdir())])
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    build=ROOT/'.local/adsorption/column_transport_domain';metadata=build_case(HERE/'openlb_domain',build,base/'build')
    prior=None
    if args.reuse_native:
        prior_path=args.reuse_native.resolve()/'summary.json'
        prior=json.loads(prior_path.read_text())
        for name,digest in prior['artifacts_sha256'].items():
            if sha(prior_path.parent/name)!=digest:raise ValueError('Saved native diagnostic artifact changed: '+name)
        for name,digest in prior['source_sha256'].items():
            if sha(prior_path.parent/'source_snapshot'/name)!=digest:
                raise ValueError('Archived diagnostic source changed: '+name)
        if not prior['checks_passed'] or prior['reference_field_sha256']!=sha(reference/'field.csv'):
            raise ValueError('Saved native diagnostics/reference are not identical qualified inputs.')
        for name,digest in metadata['source_sha256'].items():
            if prior['build']['source_sha256'].get(name)!=digest:
                raise ValueError('Native application changed; execute new diagnostic simulations.')
        shutil.copytree(prior_path.parent/'source_snapshot',base/'prior_source_snapshot')
    records={}
    for name,cells,field,perturb in [('cold',1,'none',0.),('domain',2,str(reference/'field.csv'),.05)]:
        dest=base/name
        command=[str(build/'capillary3d'),'32',repr(r['gradient_Pa_m']),repr(r['tau_even']),'.0008',field,'1.',str(cells),repr(perturb)]
        request=dict(command=command,binary_sha256=sha(build/'capillary3d'),
            initial_field_sha256=None if field=='none' else sha(Path(field)))
        if prior:
            source=args.reuse_native.resolve()/name
            saved=json.loads((source/'request.json').read_text())
            actual=list(saved['command']);expected=list(command);actual[5]=expected[5]
            if actual!=expected or any(saved[k]!=request[k] for k in ('binary_sha256','initial_field_sha256')):
                raise ValueError('Saved native argv/binary/initial field changed.')
            shutil.copytree(source,dest);returncode=0
            print('Reanalyze identical saved native diagnostic',name,flush=True)
        else:
            dest.mkdir();write_json(dest/'request.json',request)
            print('Native hydraulic diagnostic start',name,flush=True)
            with (dest/'stdout.log').open('w') as stream:
                result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=1800)
            returncode=result.returncode
        raw=json.loads((dest/'raw.json').read_text())
        a=np.genfromtxt(dest/'field.csv',delimiter=',',names=True)
        n=raw['resolution'];u=float(np.sum(a['ux_m_s'])/(cells*n**3))
        flux=np.array([np.sum(a['density_kg_m3'][a['i']==i]*a['ux_m_s'][a['i']==i])*raw['dx_m']**2 for i in range(cells*n)])
        hist=np.genfromtxt(dest/'history.csv',delimiter=',',names=True)
        stationarity={key:float(np.std(hist[key][-6:],ddof=1)/abs(np.mean(hist[key][-6:]))) for key in ('superficial_velocity_m_s','mean_lattice_energy')}
        delta=abs(u/old['superficial_velocity_m_s']-1)
        checks=dict(converged=bool(returncode==0 and raw['converged']),mean_velocity=delta<=1e-5,
            domain_periodicity=raw['periodicity_relative_error']<=1e-5,
            mass=float(np.ptp(flux)/abs(np.mean(flux)))<=1e-6,
            stationarity=max(stationarity.values())<=1.0001e-7,
            Mach=raw['max_Mach']<.05,density=raw['max_density_departure']<1e-4,
            geometry=len(a)==raw['fluid_nodes']==cells*old['raw']['fluid_nodes'])
        records[name]=dict(raw=raw,checks=checks,superficial_velocity_m_s=u,
            mean_velocity_relative_difference=delta,plane_mass_relative_spread=float(np.ptp(flux)/abs(np.mean(flux))),stationarity=stationarity)
        print(name,checks,flush=True)
    s=dict(schema='pore2process.hydraulic-domain-check.v1',run_id=rid,source_sha256=sources,
        source_reference=str(reference.relative_to(ROOT)),reference_field_sha256=sha(reference/'field.csv'),
        cases=records,checks_passed=all(all(v['checks'].values()) for v in records.values()),
        native_runs=0 if prior else 2,reused_native_runs=2 if prior else 0,
        prior_native_summary=None if prior is None else str(prior_path.relative_to(ROOT)),
        prior_native_summary_sha256=None if prior is None else sha(prior_path),
        reuse_reason=None if prior is None else 'Portable wrapper reference selection; native binary, physical argv and input bytes are identical.',
        new_DWSIM_runs=0,build=metadata,vendor=verify_vendor())
    s['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',s)
    write_json(OUTPUT/'hydraulic_diagnostics/latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Hydraulic diagnostics saved',base,s['checks_passed'],flush=True)


if __name__=='__main__':main()
