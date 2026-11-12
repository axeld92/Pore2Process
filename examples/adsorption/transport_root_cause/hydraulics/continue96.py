"""Continue the interrupted N96 pressure root without repeating its cold solve."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

from pore2process_demo.io import sha,write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest,verified_source_path
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from .run import HERE,OUTPUT,analyze,native_checks,source_point


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interrupted',type=Path,required=True)
    args=parser.parse_args();old=args.interrupted.resolve();plan=json.loads((old/'plan.json').read_text())
    if plan['resolution']!=96 or (old/'summary.json').exists():raise ValueError('Expected interrupted N96 parent')
    for path,digest in plan['source_sha256'].items():
        if sha(verified_source_path(path,digest))!=digest:raise ValueError('Original source changed: '+path)
    for role,key in [('hydraulic','hydraulic_summary_sha256'),('native_crosscheck','native_crosscheck_sha256'),('previous_refinement','previous_refinement_sha256')]:
        path=ROOT/plan[role+'_summary']
        if sha(path)!=plan[key]:raise ValueError('Original dependency changed: '+role)
        validate_record(path)
    history=json.loads((old/'grid96_progress.json').read_text())['iterations']
    if len(history)!=1:raise ValueError('Expected one fully completed pressure trial')
    cold=old/history[0]['directory'];partial=cold.parent/'iteration2'
    if (partial/'raw.json').exists():raise ValueError('Second trial already completed; do not repeat it')
    a=analyze(cold)
    if a!=history[0]['analysis']:raise ValueError('Cold result does not replay')
    request=json.loads((partial/'request.json').read_text());warm=cold/'field.csv'
    binary_source=old/'binary/capillary3d'
    if sha(warm)!=request['initial_field_sha256'] or sha(binary_source)!=plan['binary_sha256']:
        raise ValueError('Preserved field/executable changed')
    carrier=plan['inputs'];target=carrier['native_target_velocity_m_s'];u=a['superficial_velocity_m_s']
    cold_checks=native_checks(a,target)
    if not a['raw']['converged'] or not all(v for k,v in cold_checks.items() if k!='flow'):
        raise ValueError('Cold parent failed its non-flow screens')
    # The incomplete trial provides a pressure guess only, never an accepted
    # velocity field, steady-flow observation or replacement for convergence.
    last=(partial/'history.csv').read_text().splitlines()[-1].split(',')
    g0=a['raw']['gradient_Pa_m'];g1=float(request['command'][2]);u1=float(last[1])
    gradient=g1+(target-u1)*(g1-g0)/(u1-u)
    if not .5*g1<gradient<2*g1:raise ValueError('Invalid continuation pressure guess')
    interrupted=dict(plan,status='interrupted_by_host_session_update',checks_passed=False,source_point=None,
        new_native_flow_runs=1,native=dict(iterations=history),native_started=2,native_completed=1,
        native_interrupted=1,interrupted_returncode=None,
        interruption_evidence='Update ended prior session; process check found no matching OpenLB/controller. Second trial has no final field or raw result.',
        last_partial_step=int(last[0]),last_partial_superficial_velocity_m_s=u1)
    parent_path=old/'interruption.json'
    if parent_path.exists():
        if validate_record(parent_path)['last_partial_step']!=int(last[0]):raise ValueError('Preserved interruption changed')
    else:
        interrupted['artifacts_sha256']={str(p.relative_to(old)):sha(p) for p in old.rglob('*') if p.is_file()}
        write_json(parent_path,interrupted)
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'continuations'/rid;base.mkdir(parents=True)
    sources=source_manifest([Path(__file__),HERE/'run.py',HERE/'protocol.md',HERE.parent/'resolution96_protocol.md',*list((HERE/'openlb').iterdir())])
    for path in sources:
        saved=base/'source_snapshot'/path;saved.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/path,saved)
    binary=base/'binary/capillary3d';binary.parent.mkdir();shutil.copy2(binary_source,binary)
    record=dict(plan,run_id=rid,source_sha256=sources,parent_summary=str(parent_path.relative_to(ROOT)),parent_summary_sha256=sha(parent_path),
        continuation=dict(verified_cold_field_sha256=sha(warm),partial_pressure_guess_only=True,first_gradient_Pa_m=gradient,
            original_cold_solves_reused=1,original_partial_solves_reused_as_fields=0,unchanged_acceptance=True))
    write_json(base/'plan.json',record)
    env=dict(os.environ,OMP_NUM_THREADS=str(plan['OMP_NUM_THREADS']),OMP_DYNAMIC='FALSE')
    done=[];previous=(g0,u);warm_scale=target/u
    for i in range(1,7):
        dest=base/f'native_flow/d800_grid96/iteration{i}';dest.mkdir(parents=True)
        command=[str(binary),'96',repr(gradient),request['command'][3],request['command'][4],str(warm),repr(warm_scale)]
        write_json(dest/'request.json',dict(command=command,binary_sha256=sha(binary),initial_field_sha256=sha(warm),target_superficial_velocity_m_s=target))
        print('N96 continuation starts',base,i,'gradient',gradient,flush=True);start=perf_counter()
        with (dest/'stdout.log').open('w') as log:
            native=subprocess.run(command,cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=10800)
        write_json(dest/'execution.json',dict(returncode=native.returncode,elapsed_s=perf_counter()-start))
        if native.returncode:raise RuntimeError('Native pressure continuation failed')
        a=analyze(dest);u=a['superficial_velocity_m_s'];checks=native_checks(a,target)
        item=dict(directory=str(dest.relative_to(base)),analysis=a,flow_relative_error=abs(u/target-1),checks=checks,initial_field_sha256=sha(warm))
        done.append(item);write_json(base/'grid96_progress.json',dict(resolution=96,iterations=done))
        if all(checks.values()):break
        if any(not v for k,v in checks.items() if k!='flow'):raise ValueError('Pressure trial failed non-flow screens')
        old_g,old_u=previous;new_g=gradient*target/u
        if u!=old_u:
            secant=gradient+(target-u)*(gradient-old_g)/(u-old_u)
            if .5*gradient<secant<2*gradient:new_g=secant
        previous=(gradient,u);gradient=new_g;warm=dest/'field.csv';warm_scale=target/u
    else:raise RuntimeError('Fixed flow not reached within original pressure trial budget')
    point=source_point(dest,carrier);write_json(base/'source_point.json',point);r=a['raw']
    record.update(native=dict(accepted=item,iterations=done,checks=checks),source_point=point,new_native_flow_runs=len(done),
        geometry=dict(resolution=96,fluid_nodes=r['fluid_nodes'],side_m=r['unit_cell_side_m'],dx_m=r['dx_m'],
            domain_volume_m3=r['unit_cell_side_m']**3,physical_liquid_volume_m3=r['fluid_nodes']*r['dx_m']**3,
            voxel_porosity=point['voxel_porosity'],particle_diameter_m=r['particle_diameter_m'],
            coordinate_system='Cartesian metres, (index+0.5)*dx; zero-based periodic i,j,k.'),
        checks_passed=True,vendor=verify_vendor())
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/'grid96_omp_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('N96 qualified',base,point,flush=True)


if __name__=='__main__':main()
