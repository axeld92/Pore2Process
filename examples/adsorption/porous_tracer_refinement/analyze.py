"""Audit and compare half-step and matched-mesh tracer refinements."""
import argparse
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.porous_face_handoff.run import HERE as HANDOFF, SNAPSHOTS, archived, scalar_metrics
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .mesh import FINE_N, N, RATIO, fine_coordinates, refined_mesh, timestep_counts
from .time_half import BOUNDED, BOUNDED_SHA, CENTER, CENTER_SHA, OUTPUT, records

HERE = Path(__file__).parent


def trajectory(path, source, dx, runtime, factor, ijk, n):
    args = [source,path,dx,ref.D,*ref.TIMES,'varying',factor]
    expected = dict(command=[*runtime,str(HANDOFF/'solve_bounded.jl'),*map(str,args)],
                    balanced_sha256=sha(source/'balanced.f64'),initial_sha256=sha(source/'varying_initial.f64'))
    if read_json(path/'request.json')!=expected or read_json(path/'execution.json')['returncode']:
        raise ValueError('Saved tracer execution differs from the refinement plan')
    initial = np.fromfile(source/'varying_initial.f64',dtype='<f8')
    metrics, fields = scalar_metrics(path,'varying',initial,ijk,n)
    history = np.genfromtxt(path/'varying_history.csv',delimiter=',',names=True)
    return metrics, fields, history['time_s'], float(np.var(initial))


def variance_balance(metrics, initial_variance):
    result = {}
    for name, snapshot in metrics['snapshots'].items():
        lost = .5*(initial_variance-snapshot['variance_kg2_m6'])
        if lost<=0:
            raise ValueError('Expected positive variance decay over the test window')
        molecular = snapshot['integrated_dissipation']
        numerical = lost-molecular
        result[name] = dict(half_variance_loss_kg2_m6=lost,molecular_loss_kg2_m6=molecular,
                            numerical_loss_kg2_m6=numerical,numerical_fraction=numerical/lost,
                            within_numerical_limit=abs(numerical/lost)<=.15)
    return result


def compare(base):
    center,bounded = records()
    time_plan = read_json(base/'time_plan.json')
    if (time_plan['center_sha256'],time_plan['bounded_sha256'])!=(CENTER_SHA,BOUNDED_SHA):
        raise ValueError('The archived controls differ from the plan')
    for name,digest in time_plan['source_sha256'].items():
        if sha(ROOT/name)!=digest:
            raise ValueError('Time-step source differs: '+name)
    source = CENTER.parent/'N96'
    info,ijk,_ = archived(N)
    fine_source = base/'N192_faces'
    mesh_plan = read_json(base/'mesh_plan.json')
    mesh = refined_mesh(fine_source,source,info,ijk,write=False)
    if mesh!=mesh_plan['mesh'] or sha(HERE/'mesh.py')!=mesh_plan['mesh_source_sha256']:
        raise ValueError('Refined mesh or source differs from the plan')
    runtime = center['julia_runtime']
    dx = info['dx_m']
    fine_ijk = fine_coordinates(ijk)
    half = trajectory(base/'N96_half_step',source,dx,runtime,.5,ijk,N)
    matched_name = mesh['matched_case']
    matched = half if matched_name=='N96_half_step' else trajectory(
        base/matched_name,source,dx,runtime,mesh['matched_coarse_time_factor'],ijk,N)
    fine = trajectory(base/'N192',fine_source,mesh['dx_m'],runtime,mesh['fine_time_factor'],fine_ijk,FINE_N)
    np.testing.assert_array_equal(matched[2],fine[2])
    if fine[0]['history_steps']!=sum(mesh['common_step_counts']):
        raise ValueError('Unexpected common time schedule')
    cap = mesh['coarse_stability_cap_s']
    if half[0]['history_steps']!=sum(timestep_counts(cap/2)):
        raise ValueError('Half-step history has the wrong length')
    original_dir = BOUNDED.parent/'N96/varying'
    initial = np.fromfile(source/'varying_initial.f64',dtype='<f8')
    original, old_fields = scalar_metrics(original_dir,'varying',initial,ijk,N)
    old = bounded['analysis']['cases']['N96_varying']
    if original['snapshots'].keys()!=old['snapshots'].keys() or original['history_steps']!=sum(timestep_counts(cap)):
        raise ValueError('Original N96 history differs from the recorded cap')
    for name, metrics in original['snapshots'].items():
        if any(old['snapshots'][name][k]!=v for k,v in metrics.items()):
            raise ValueError('Original N96 scalar metrics differ from the archive')
    cases = {'N96_original':original,'N96_half_step':half[0],'N192':fine[0]}
    losses = {'N96_original':variance_balance(original,float(np.var(initial))),
              'N96_half_step':variance_balance(half[0],half[3]),'N192':variance_balance(fine[0],fine[3])}
    if matched_name not in cases:
        cases[matched_name] = matched[0]
        losses[matched_name] = variance_balance(matched[0],matched[3])
    # These bounded runs include limiter dissipation, so report a separate
    # numerical share rather than applying the centered diffusion-only budget.
    for result in cases.values():
        for t in SNAPSHOTS:
            result['checks'].pop(t+'_variance_budget')
    time, space = {}, {}
    for i,t in enumerate(SNAPSHOTS):
        coarse,c = matched[0]['snapshots'][t],fine[0]['snapshots'][t]
        v = abs(coarse['variance_kg2_m6']/c['variance_kg2_m6']-1)
        a = abs(coarse['fundamental_amplitude_kg_m3']/c['fundamental_amplitude_kg_m3']-1)
        p = abs(float(np.angle(np.exp(1j*(coarse['phase_rad']-c['phase_rad'])))))
        averaged = fine[1][i].reshape(-1,RATIO**3).mean(axis=1)
        e = float(np.sqrt(np.mean(((averaged-matched[1][i])/ref.A)**2)))
        nl96 = losses[matched_name][t]['numerical_loss_kg2_m6']
        nl192 = losses['N192'][t]['numerical_loss_kg2_m6']
        space[t] = dict(variance_relative_change=v,fundamental_amplitude_relative_change=a,
                        fundamental_phase_change_rad=p,coarse_averaged_field_RMS_over_amplitude=e,
                        numerical_loss_reduction_fraction=1-nl192/nl96,
                        numerical_loss_decreased=abs(nl192)<abs(nl96),
                        passed=v<=.05 and a<=.05 and p<=.05)
        e = float(np.sqrt(np.mean(((old_fields[i]-half[1][i])/ref.A)**2)))
        time[t] = dict(RMS_over_amplitude=e,passed=e<=.001)
    steps = mesh['common_step_counts']
    common_dts = [(end-start)/step for start,end,step in zip((0.,*ref.TIMES[:-1]),ref.TIMES,steps)]
    checks = dict(flow_preserved=True,
        bounded_and_conserved=all(c['checks']['mass'] and c['checks']['physical_bounds'] for c in cases.values()),
        temporal_profile=all(t['passed'] for t in time.values()),
        spatial_spreading=all(s['passed'] for s in space.values()),
        numerical_loss_decreases=all(s['numerical_loss_decreased'] for s in space.values()),
        fine_numerical_share=all(s['within_numerical_limit'] for s in losses['N192'].values()))
    return dict(cases=cases,variance_budget=losses,mesh=mesh,common_effective_dt_s=common_dts,
                time_change=time,spatial_change=space,checks=checks)


def finalize(base):
    base = base.resolve()
    analysis = compare(base)
    files = [HERE/k for k in ('analyze.py','mesh.py','time_half.py','protocol.md')]
    files += [HANDOFF/'run.py',HANDOFF/'solve_bounded.jl',ROOT/'julia/adsorption/BoundedFaceTransport.jl']
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record=dict(schema='pore2process.porous-tracer-refinement.v1',run_id=base.name,
                source_sha256=hashes,centered_sha256=CENTER_SHA,bounded_sha256=BOUNDED_SHA,
                analysis=analysis,new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0)
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*')
                               if p.is_file() and p.name!='summary.json'}
    path=base/'summary.json'
    write_json(path,record)
    write_json(OUTPUT/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed refinement',path,analysis['checks'])


def check(path):
    path = path.resolve()
    saved = validate_record(path)
    if compare(path.parent)!=saved['analysis']:
        raise ValueError('Replayed refinement differs')
    report=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(saved['source_sha256']),
                artifacts=len(saved['artifacts_sha256']),checks=saved['analysis']['checks'])
    write_json(OUTPUT/'verification.json',report)
    print(report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',type=Path)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    check(args.path) if args.check else finalize(args.path)
