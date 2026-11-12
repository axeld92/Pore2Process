"""Run the existing bounded Julia tracer at half its original N96 time cap."""
from datetime import datetime, timezone
from pathlib import Path
import subprocess
from time import perf_counter
from uuid import uuid4

from examples.adsorption.porous_face_handoff.run import HERE as HANDOFF, OUTPUT as PARENT_OUTPUT
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT

OUTPUT = ROOT/'results/adsorption/porous_tracer_refinement'
CENTER = PARENT_OUTPUT/'runs/20261111T154427Z-6e178a5c/summary.json'
BOUNDED = PARENT_OUTPUT/'bounded/20261111T154712Z-4a96ae28/summary.json'
CENTER_SHA = 'd3452087e2c39da82440dac6d8aab47477220f41a15ca72fbb813f505bf028cd'
BOUNDED_SHA = '7968ae504da137c9b5be23d7c8e2f09615c8ac856debf667848c0d0c4055faee'


def records():
    if sha(CENTER)!=CENTER_SHA or sha(BOUNDED)!=BOUNDED_SHA:
        raise ValueError('The archived N96 controls changed')
    center, bounded = validate_record(CENTER), validate_record(BOUNDED)
    if bounded['centered_summary_sha256'] != CENTER_SHA or not center['analysis']['handoff_qualified']:
        raise ValueError('The bounded study did not use the qualified saved faces')
    return center, bounded


def run_solver(directory, source, dx, runtime, time_factor):
    directory.mkdir(parents=True)
    args = [source,directory,dx,ref.D,*ref.TIMES,'varying',time_factor]
    cmd = [*runtime,str(HANDOFF/'solve_bounded.jl'),*map(str,args)]
    write_json(directory/'request.json',dict(command=cmd,balanced_sha256=sha(source/'balanced.f64'),
                                              initial_sha256=sha(source/'varying_initial.f64')))
    print('Julia tracer start',directory,flush=True)
    start = perf_counter()
    with (directory/'stdout.log').open('w') as log:
        process = subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    write_json(directory/'execution.json',dict(returncode=process.returncode,wall_s=perf_counter()-start))
    if process.returncode:
        raise RuntimeError(f'Julia tracer failed. See {directory / "stdout.log"}')
    print('Julia tracer complete',directory,round(perf_counter()-start,1),'s',flush=True)


def run():
    center, _ = records()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'runs'/stamp
    base.mkdir(parents=True)
    source = CENTER.parent/'N96'
    dx = center['analysis']['cases']['96']['source']['dx_m']
    write_json(base/'time_plan.json',dict(center_sha256=CENTER_SHA,bounded_sha256=BOUNDED_SHA,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (
            Path(__file__),Path(__file__).with_name('protocol.md'),HANDOFF/'solve_bounded.jl',
            ROOT/'julia/adsorption/BoundedFaceTransport.jl')},
        new_OpenLB_runs=0,new_bed_runs=0,new_DWSIM_runs=0))
    print('Refinement study',base,flush=True)
    run_solver(base/'N96_half_step',source,dx,center['julia_runtime'],.5)


if __name__=='__main__':
    run()
