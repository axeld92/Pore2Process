"""Independently solve and compare a translated same-material FCC problem."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
from uuid import uuid4
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from examples.adsorption.bed import analyze as analyze_water
from examples.adsorption.transport_accuracy.run import validate_record
from ..hydraulics.run import refine,source_point
from ..study import DT,END,execute,analyze as analyze_scalar_base
from ..flux_audit.analysis import summarize,summarize_uniform_stencil
from .geometry import sphere_distance,location

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/transport_root_cause/grid_phase'
BASELINE=ROOT/'results/adsorption/transport_root_cause/runs/20261110T205424Z-ed34b563/summary.json'


def analyze_scalar(path,phase_x_m):
    path=Path(path);r=analyze_scalar_base(path);raw=r['raw']
    f=np.genfromtxt(path/'peak.csv',delimiter=',',names=True);fluid=f['material']!=2
    xyz=np.column_stack([f[k+'_m'] for k in ('x','y','z')])
    gaps,_,_=sphere_distance(xyz,side_m=raw['side_m'],phase_x_m=phase_x_m)
    if not np.array_equal(fluid,gaps>0):raise ValueError('Exported mask differs from analytical translated spheres.')
    ids=np.flatnonzero(fluid)
    extrema={name:location(f[index],side_m=raw['side_m'],dx_m=raw['dx_m'],phase_x_m=phase_x_m)
        for name,index in [('minimum',ids[np.argmin(f['concentration_kg_m3'][fluid])]),
                          ('maximum',ids[np.argmax(f['concentration_kg_m3'][fluid])]) ]}
    audits={p.stem:summarize(p,dx_m=raw['dx_m'],dt_s=raw['dt_s']) for p in sorted(path.glob('audit_*.csv'))}
    return dict(result=r,phase_x_m=phase_x_m,fluid_nodes=int(np.sum(fluid)),
        voxel_porosity=float(np.mean(fluid)),analytical_mask_matches=True,
        minimum_fluid_gap_m=float(np.min(gaps[fluid])),extrema=extrema,
        exact_population_audits=audits,uniform_stencil=summarize_uniform_stencil(path/'uniform_stencil.csv'))


def prefix_comparison(full,prefix):
    full=Path(full);prefix=Path(prefix)
    a=np.genfromtxt(full/'peak.csv',delimiter=',',names=True);b=np.genfromtxt(prefix/'peak.csv',delimiter=',',names=True)
    field=max(float(np.max(abs(a[k]-b[k]))) for k in a.dtype.names)
    a=np.genfromtxt(full/'inventory.csv',delimiter=',',names=True);b=np.genfromtxt(prefix/'inventory.csv',delimiter=',',names=True)
    if len(b)>len(a):raise ValueError('Prefix extends beyond full run.')
    inventory=max(float(np.max(abs(a[k][:len(b)]-b[k]))) for k in a.dtype.names)
    return dict(maximum_peak_field_difference=field,maximum_inventory_prefix_difference=inventory,
                passed=field<=1e-12 and inventory<=1e-12)


def compare(unshifted,shifted,old_water,new_water):
    a,b=unshifted['result'],shifted['result']
    if abs(a['raw']['dt_s']/b['raw']['dt_s']-1)>1e-12:raise ValueError('Grid-phase dt differs.')
    if abs(a['actual_duration_s']-b['actual_duration_s'])>1e-12:raise ValueError('Grid-phase duration differs.')
    old,new=a['peak']['maximum_excursion_kg_m3'],b['peak']['maximum_excursion_kg_m3']
    return dict(excursion_reduction_fraction=None if old==0 else 1-new/old,
        unshifted_voxel_porosity=unshifted['voxel_porosity'],shifted_voxel_porosity=shifted['voxel_porosity'],
        voxel_porosity_difference=shifted['voxel_porosity']-unshifted['voxel_porosity'],
        pressure_gradient_relative_change=new_water['raw']['gradient_Pa_m']/old_water['raw']['gradient_Pa_m']-1,
        apparent_operating_point_permeability_relative_change=new_water['permeability_m2']/old_water['permeability_m2']-1,
        superficial_velocity_relative_change=new_water['superficial_velocity_m_s']/old_water['superficial_velocity_m_s']-1,
        final_liquid_SD_relative_change=b['periodic_wave_final']['sd_C_kg_m3']/a['periodic_wave_final']['sd_C_kg_m3']-1,
        final_first_mode_amplitude_relative_change=np.hypot(b['periodic_wave_final']['sin_coefficient_kg_m3'],b['periodic_wave_final']['cos_coefficient_kg_m3'])/np.hypot(a['periodic_wave_final']['sin_coefficient_kg_m3'],a['periodic_wave_final']['cos_coefficient_kg_m3'])-1,
        interpretation='Alignment sensitivity of translated continuum problem; discrete voxel geometry changes. Not material improvement or a universal fix.')


def replay(path):
    p=Path(path);s=validate_record(p);base=p.parent
    baseline=ROOT/s['baseline_summary'];old=validate_record(baseline)
    if sha(baseline)!=s['baseline_summary_sha256']:raise ValueError('Grid-phase baseline changed.')
    hydraulic=ROOT/s['hydraulic_summary'];validate_record(hydraulic)
    if sha(hydraulic)!=s['hydraulic_summary_sha256']:raise ValueError('Grid-phase hydraulic parent changed.')
    old_dest=baseline.parent/'native/n32_sine'
    if analyze_scalar(old_dest,0.)!=s['unshifted']:raise ValueError('Unshifted analysis differs.')
    for key,value in s['shifted'].items():
        if analyze_scalar(base/'native'/key,s['phase_x_m'])!=value:raise ValueError('Shifted analysis differs: '+key)
    if prefix_comparison(base/'native/shifted_sine',base/'native/shifted_peak_audit')!=s['prefix_reproduction']:raise ValueError('Peak prefix differs.')
    water=base/s['native_flow']['accepted']['directory'];new_water=analyze_water(water)
    if source_point(water,s['inputs'])!=s['source_point']:raise ValueError('Shifted water point differs.')
    raw=new_water['raw']
    if summarize(water/'water_population_audit.csv',dx_m=raw['dx_m'],dt_s=raw['dt_s'],density_scale=1000.)!=s['exact_water_population_audit']:
        raise ValueError('Shifted exact water audit differs.')
    old_point=next(item for item in old['specifications'] if item['id']=='n32_sine');old_water=analyze_water(Path(old_point['field']).parent)
    if compare(s['unshifted'],s['shifted']['shifted_sine'],old_water,new_water)!=s['comparison']:raise ValueError('Grid-phase comparison differs.')
    return s


def check():
    p=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];s=replay(p)
    record=dict(summary=str(p.relative_to(ROOT)),summary_sha256=sha(p),run_id=s['run_id'],
                source_files_checked=len(s['source_sha256']),artifacts_checked=len(s['artifacts_sha256']),
                actual_native_flow_runs=s['new_native_flow_runs'],actual_native_scalar_runs=s['new_native_scalar_runs'],
                new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,vendor=verify_vendor(),
                comparison=s['comparison'],prefix_reproduction=s['prefix_reproduction'],
                unshifted=s['unshifted'],shifted=s['shifted'])
    write_json(OUTPUT/'verification.json',record);write_json(HERE/'reference.json',record)
    print(json.dumps(record,indent=2));return record


def run(baseline=BASELINE):
    baseline=Path(baseline).resolve();old=validate_record(baseline)
    old_point=next(item for item in old['specifications'] if item['id']=='n32_sine')
    hydraulic=ROOT/old['hydraulic_summary'];hp=validate_record(hydraulic)
    if sha(hydraulic)!=old['hydraulic_summary_sha256']:raise ValueError('Old hydraulic parent changed.')
    carrier=json.loads((hydraulic.parent/'inputs.json').read_text());old_water=analyze_water(Path(old_point['field']).parent)
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8];base=OUTPUT/'runs'/rid;base.mkdir(parents=True)
    files=[HERE/'run.py',HERE/'geometry.py',HERE/'protocol.md',*list((HERE/'openlb_flow').iterdir()),*list((HERE/'openlb_scalar').iterdir()),
           HERE.parent/'study.py',HERE.parent/'hydraulics/run.py',HERE.parent/'flux_audit/analysis.py',
           ROOT/'examples/adsorption/bed.py',ROOT/'examples/adsorption/transport_accuracy/profile.py',ROOT/'examples/adsorption/transport_accuracy/run.py',ROOT/'examples/adsorption/column_transport/run.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    phase=old_point['side_m']/64
    record=dict(schema='pore2process.fcc-grid-phase-diagnostic.v1',run_id=rid,source_sha256=sources,
        baseline_summary=str(baseline.relative_to(ROOT)),baseline_summary_sha256=sha(baseline),
        hydraulic_summary=str(hydraulic.relative_to(ROOT)),hydraulic_summary_sha256=sha(hydraulic),
        phase_x_m=phase,inputs=carrier,scalar_dt_s=DT,scalar_requested_end_s=END,
        new_DWSIM_runs=0,new_adsorption_predictions=0,service_prediction_qualified=False)
    write_json(base/'plan.json',record);print('Grid-phase run',base,flush=True)
    build_flow=build_case(HERE/'openlb_flow',ROOT/'.local/adsorption/grid_phase_flow',base/'build_flow')
    binaries=base/'binary';binaries.mkdir();flow_binary=binaries/'flow';shutil.copy2(ROOT/'.local/adsorption/grid_phase_flow/capillary3d',flow_binary)
    group=refine(base,flow_binary,32,carrier['native_target_velocity_m_s'],old_water['raw']['gradient_Pa_m'])
    water=base/group['accepted']['directory'];point=source_point(water,carrier);write_json(base/'source_point.json',point)
    build_scalar=build_case(HERE/'openlb_scalar',ROOT/'.local/adsorption/grid_phase_scalar',base/'build_scalar')
    scalar_binary=binaries/'scalar';shutil.copy2(ROOT/'.local/adsorption/grid_phase_scalar/capillary3d',scalar_binary)
    spec=dict(point,id='shifted_sine',n=32,phase_x_m=phase,amplitude=.25,audit_step=-1,end_s=END)
    execute(base,spec,scalar_binary);full=analyze_scalar(base/'native/shifted_sine',phase)
    peak_step=full['result']['peak']['step'];prefix=dict(spec,id='shifted_peak_audit',audit_step=peak_step,
        end_s=(peak_step-.25)*full['result']['raw']['dt_s'])
    execute(base,prefix,scalar_binary);reproduction=prefix_comparison(base/'native/shifted_sine',base/'native/shifted_peak_audit')
    if not reproduction['passed']:raise ValueError('Read-only peak prefix does not reproduce full run.')
    unshifted=analyze_scalar(baseline.parent/'native/n32_sine',0.)
    flow_raw=group['accepted']['analysis']['raw']
    record.update(build_flow=build_flow,build_scalar=build_scalar,native_flow=group,source_point=point,
        new_native_flow_runs=len(group['iterations']),new_native_scalar_runs=2,
        specifications=[spec,prefix],unshifted=unshifted,
        shifted={name:analyze_scalar(base/'native'/name,phase) for name in ('shifted_sine','shifted_peak_audit')},
        prefix_reproduction=reproduction,comparison=compare(unshifted,full,old_water,analyze_water(water)),
        exact_water_population_audit=summarize(water/'water_population_audit.csv',dx_m=flow_raw['dx_m'],dt_s=flow_raw['dt_s'],density_scale=1000.),
        vendor=verify_vendor())
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record);write_json(OUTPUT/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Grid-phase complete',base,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');p.add_argument('--baseline',type=Path,default=BASELINE);args=p.parse_args()
    check() if args.check else run(args.baseline)

if __name__=='__main__':main()
