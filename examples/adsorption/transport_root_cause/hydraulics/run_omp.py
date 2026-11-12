"""Use an explicitly cross-checked stock OpenMP binary for fixed-flow grids."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
from uuid import uuid4

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from .run import HERE,OUTPUT,latest_parent,refine,source_point


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--crosscheck',type=Path,required=True)
    parser.add_argument('--grid',type=int,choices=[64,96],default=64)
    parser.add_argument('--previous-refinement',type=Path)
    args=parser.parse_args()
    cross_path=args.crosscheck.resolve()/'summary.json'
    cross=validate_record(cross_path)
    if not cross['checks_passed'] or cross['measured_speedup']<=1:
        raise ValueError('A successful, faster native serial/OpenMP cross-check is required.')
    hydraulic=latest_parent().resolve()
    parent=validate_record(hydraulic/'summary.json')
    if cross['hydraulic_summary_sha256']!=sha(hydraulic/'summary.json'):
        raise ValueError('OpenMP cross-check used another physical operating point.')
    carrier=parent['inputs']
    previous=None
    gradient=parent['cases']['d800']['native']['48']['accepted']['analysis']['raw']['gradient_Pa_m']
    if args.previous_refinement:
        previous=validate_record(args.previous_refinement.resolve()/'summary.json')
        if previous['hydraulic_summary_sha256']!=sha(hydraulic/'summary.json'):
            raise ValueError('Previous refinement used another operating point.')
        gradient=previous['native']['accepted']['analysis']['raw']['gradient_Pa_m']
    binary_source=Path(cross['binary'])
    if sha(binary_source)!=cross['binary_sha256']:
        raise ValueError('Cross-checked native OpenMP executable changed.')
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'omp_runs'/rid;base.mkdir(parents=True)
    files=[Path(__file__),HERE/'run.py',HERE/'omp_check.py',HERE/'protocol.md',
        *list((HERE/'openlb').iterdir()),
        ROOT/'examples/adsorption/transport_root_cause/flux_audit/population_audit.h',
        ROOT/'examples/adsorption/bed.py',ROOT/'examples/adsorption/transport_accuracy/run.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    binary=base/'binary/capillary3d';binary.parent.mkdir();shutil.copy2(binary_source,binary)
    os.environ['OMP_NUM_THREADS']=str(cross['threads'])
    os.environ['OMP_DYNAMIC']='FALSE'
    write_json(base/'inputs.json',carrier)
    write_json(base/'controller.json',dict(pid=os.getpid(),scope='This fixed-flow native run only.'))
    record=dict(schema='pore2process.transport-root-hydraulic-refinement.v1',run_id=rid,
        resolution=args.grid,source_sha256=sources,
        hydraulic_summary=str((hydraulic/'summary.json').relative_to(ROOT)),
        hydraulic_summary_sha256=sha(hydraulic/'summary.json'),inputs=carrier,
        parallel_mode='stock OpenLB CPU_SISD OpenMP',OMP_NUM_THREADS=cross['threads'],OMP_DYNAMIC=False,
        native_crosscheck_summary=str(cross_path.relative_to(ROOT)),native_crosscheck_sha256=sha(cross_path),
        previous_refinement_summary=None if previous is None else str((args.previous_refinement.resolve()/'summary.json').relative_to(ROOT)),
        previous_refinement_sha256=None if previous is None else sha(args.previous_refinement.resolve()/'summary.json'),
        build=cross['build'],binary_sha256=sha(binary),
        new_DWSIM_runs=0,new_adsorption_predictions=0,
        geometry_role='Unchanged stock synthetic FCC sphere indicators at diameter0.8mm nominal porosity0.6.')
    write_json(base/'plan.json',record)
    print('OpenMP hydraulic refinement',base,flush=True)
    group=refine(base,binary,args.grid,carrier['native_target_velocity_m_s'],gradient)
    accepted=base/group['accepted']['directory']
    point=source_point(accepted,carrier);write_json(base/'source_point.json',point)
    raw=group['accepted']['analysis']['raw']
    record.update(native=group,source_point=point,new_native_flow_runs=len(group['iterations']),
        geometry=dict(resolution=args.grid,fluid_nodes=raw['fluid_nodes'],
            side_m=raw['unit_cell_side_m'],dx_m=raw['dx_m'],domain_volume_m3=raw['unit_cell_side_m']**3,
            physical_liquid_volume_m3=raw['fluid_nodes']*raw['dx_m']**3,
            voxel_porosity=point['voxel_porosity'],particle_diameter_m=raw['particle_diameter_m'],
            coordinate_system='Cartesian metres, (index+0.5)*dx; zero-based periodic i,j,k.'),
        checks_passed=all(group['checks'].values()),vendor=verify_vendor())
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',record)
    write_json(OUTPUT/f'grid{args.grid}_omp_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('OpenMP hydraulic refinement accepted',base,point,flush=True)


if __name__=='__main__':main()
