"""Halve the N192 effective time step with all archived tracer inputs fixed."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.porous_face_handoff.run import HERE as HANDOFF, SNAPSHOTS, archived
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .analyze import trajectory, variance_balance
from .mesh import N, FINE_N, fine_coordinates, timestep_counts
from .time_half import OUTPUT, records, run_solver

HERE = Path(__file__).parent
PARENT = OUTPUT/'runs/20261111T160638Z-156b2680/summary.json'
PARENT_SHA = '3f4de99048825a59235c8193d2d5d357b3cf7ec0246a7bab4e6226544e30e0ef'
RESULTS = OUTPUT/'N192_time'


def inputs():
    if sha(PARENT)!=PARENT_SHA:
        raise ValueError('Archived N192 summary changed')
    parent = validate_record(PARENT)
    center,_ = records()
    mesh = parent['analysis']['mesh']
    counts = [2*n for n in mesh['common_step_counts']]
    widths = np.diff((0.,*ref.TIMES))
    # Choose an interior cap in the intersection of the intervals that
    # make ceil(window/cap) equal to each requested count.
    lower = max(float(w/n) for w,n in zip(widths,counts,strict=True))
    upper = min(float(w/(n-1)) for w,n in zip(widths,counts,strict=True))
    if not 0 < lower < upper <= mesh['fine_stability_cap_s']:
        raise ValueError('Cannot represent exact halving with a common stable cap')
    cap = (lower+upper)/2
    if timestep_counts(cap)!=counts:
        raise ValueError('Time cap does not yield exact step doubling')
    plan = dict(parent_summary=str(PARENT.relative_to(ROOT)),parent_summary_sha256=PARENT_SHA,
        old_step_counts=mesh['common_step_counts'],new_step_counts=counts,
        old_dt_s=parent['analysis']['common_effective_dt_s'],new_dt_s=[float(w/n) for w,n in zip(widths,counts,strict=True)],
        dx_m=mesh['dx_m'],dt_cap_s=cap,time_factor=cap/mesh['fine_stability_cap_s'])
    return parent,center['julia_runtime'],plan


def compare(base,parent,runtime,plan):
    _,coarse_ijk,_ = archived(N)
    ijk = fine_coordinates(coarse_ijk)
    source = PARENT.parent/'N192_faces'
    baseline, new = (trajectory(path,source,plan['dx_m'],runtime,factor,ijk,FINE_N)
        for path,factor in [(PARENT.parent/'N192',parent['analysis']['mesh']['fine_time_factor']),
                            (base/'N192_half_step',plan['time_factor'])])
    baseline_stats, baseline_fields, baseline_times, baseline_variance = baseline
    new_stats, new_fields, new_times, new_variance = new
    np.testing.assert_array_equal(new_times[::2],baseline_times)
    if new_stats['history_steps']!=sum(plan['new_step_counts']) or baseline_variance!=new_variance:
        raise ValueError('Expected exact time-step halving with the same initial field')
    budgets = dict(original=variance_balance(baseline_stats,baseline_variance),
                   half_step=variance_balance(new_stats,new_variance))
    changes = {}
    for key,old_field,new_field in zip(SNAPSHOTS,baseline_fields,new_fields,strict=True):
        old = baseline_stats['snapshots'][key]
        prior = parent['analysis']['cases']['N192']['snapshots'][key]
        if old!=prior:
            raise ValueError('Baseline metrics differ from archived N192 metrics')
        new = new_stats['snapshots'][key]
        dv = new['variance_kg2_m6']-old['variance_kg2_m6']
        dp = float(np.angle(np.exp(1j*(new['phase_rad']-old['phase_rad']))))
        rms = float(np.sqrt(np.mean(((new_field-old_field)/ref.A)**2)))
        changes[key] = dict(variance_change_kg2_m6=dv,relative_variance_change=dv/old['variance_kg2_m6'],
            fundamental_phase_change_rad=dp,
            fundamental_amplitude_relative_change=new['fundamental_amplitude_kg_m3']/old['fundamental_amplitude_kg_m3']-1,
            profile_RMS_over_initial_amplitude=rms,temporal_profile_passed=rms<=.001)
    for metrics in (baseline_stats,new_stats):
        for t in SNAPSHOTS:
            metrics['checks'].pop(t+'_variance_budget')
    return dict(original=baseline_stats,half_step=new_stats,variance_budget=budgets,changes=changes,
        prior_spatial_changes=parent['analysis']['spatial_change'],
        checks=dict(exact_half_steps=True,physical_bounds=new_stats['checks']['physical_bounds'],
                    liquid_inventory=new_stats['checks']['mass'],
                    temporal_profile=all(c['temporal_profile_passed'] for c in changes.values()),
                    numerical_share=all(b['within_numerical_limit'] for b in budgets['half_step'].values())))


def run():
    parent,runtime,plan = inputs()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = RESULTS/stamp
    base.mkdir(parents=True)
    files = [HERE/k for k in ('n192_time.py','n192_time_protocol.md','analyze.py','mesh.py','time_half.py')]
    files += [HANDOFF/'run.py',HANDOFF/'solve_bounded.jl',ROOT/'julia/adsorption/BoundedFaceTransport.jl']
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record = dict(schema='pore2process.N192-time-check.v1',run_id=stamp,plan=plan,source_sha256=hashes,
                  new_openlb_runs=0,new_meshes=0,new_bed_runs=0,new_DWSIM_runs=0)
    write_json(base/'plan.json',record)
    print('N192 time check',base,plan,flush=True)
    run_solver(base/'N192_half_step',PARENT.parent/'N192_faces',plan['dx_m'],runtime,plan['time_factor'])
    record['analysis'] = compare(base,parent,runtime,plan)
    record['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path = base/'summary.json'
    write_json(path,record)
    write_json(RESULTS/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed N192 time check',path,record['analysis']['checks'],flush=True)


def check(path):
    path = path.resolve()
    saved = validate_record(path)
    parent,runtime,plan = inputs()
    if saved['plan']!=plan or saved['analysis']!=compare(path.parent,parent,runtime,plan):
        raise ValueError('N192 time comparison differs on replay')
    result = dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(saved['source_sha256']),
                  artifacts=len(saved['artifacts_sha256']),checks=saved['analysis']['checks'])
    write_json(RESULTS/'verification.json',result)
    print(result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    check(args.check) if args.check else run()
