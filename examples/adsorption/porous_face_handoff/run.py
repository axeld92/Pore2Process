"""Transfer saved stock-TRT FCC velocities to balanced Julia tracer faces."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).parent
BASE = ROOT / 'results/adsorption'
OUTPUT = BASE / 'porous_face_handoff'
SOURCES = {
    32: ('column_transport/hydraulics/runs/20261110T154704Z-d14d13ed',
         '57ce925bd15d8bb11e0949d61f362691e9e80481ff6fb8e5fdc7d4d014e9f287'),
    64: ('transport_root_cause/hydraulics/omp_runs/20261110T210317Z-fdf047ba',
         '8692829b62c867de286717250c558e922a6cd4b4a44187fb63d3a2649f5cc2a6'),
    96: ('transport_root_cause/hydraulics/continuations/20261110T223616Z-3f4c960a',
         '25ca19605009c26b487b0f57845bef9987e4df7d9e98faccb4eb1e3ee978ce38'),
}
SNAPSHOTS = ('middle', 'final')


def archived(n):
    stem, digest = SOURCES[n]
    folder = BASE / stem
    path = folder / 'summary.json'
    if sha(path) != digest:
        raise ValueError('Saved hydraulic summary changed')
    s = read_json(path)
    inputs = s['inputs']
    conversion = (inputs['carrier_viscosity_Pa_s']/inputs['carrier_density_kg_m3']
                  /inputs['native_kinematic_viscosity_m2_s'])
    if n == 32:
        native_case = s['cases']['d800']['native'][str(n)]
        native = native_case['accepted']
        field = folder / native['directory'] / 'field.csv'
        if not all(native_case['checks'].values()):
            raise ValueError('The archived hydraulic field failed its throughput checks')
        expected = native['analysis']['superficial_velocity_m_s']*conversion
    else:
        point = s['source_point']
        field = Path(point['field'])
        if not s['checks_passed'] or point['velocity_scale'] != conversion:
            raise ValueError('Unqualified or incorrectly scaled hydraulic field')
        expected = point['actual_superficial_velocity_m_s']
    raw_file = field.with_name('raw.json')
    for p in (field, raw_file):
        if sha(p) != s['artifacts_sha256'][str(p.relative_to(folder))]:
            raise ValueError('Saved hydraulic artifact changed')
    raw = read_json(raw_file)
    if (raw['solver'], raw['collision'], raw['resolution'], raw['converged']) != ('OpenLB 1.9r0','TRT',n,True):
        raise ValueError('Expected converged stock TRT field')
    if raw['particle_diameter_m'] != .0008 or raw['unit_cell_side_m'] != ref.L:
        raise ValueError('Wrong porous geometry')
    h = raw['dx_m']
    table = np.loadtxt(field, delimiter=',', skiprows=1)
    if table.shape != (raw['fluid_nodes'],10) or not np.isfinite(table).all():
        raise ValueError('Saved field has wrong shape or non-finite values')
    ijk = table[:,:3].astype(np.int32)
    np.testing.assert_array_equal(ijk, table[:,:3])
    np.testing.assert_allclose(table[:,3:6], (ijk+.5)*h, rtol=1e-13, atol=1e-15)
    if (ijk<0).any() or (ijk>=n).any() or len(np.unique(ijk,axis=0))!=len(ijk):
        raise ValueError('Invalid Cartesian cell IDs')
    velocity = table[:,6:9]*conversion
    target = np.sum(velocity,axis=0)*h**2/n
    np.testing.assert_allclose(target[0], expected*ref.L**2, rtol=1e-12, atol=0)
    info = dict(summary=str(path.relative_to(ROOT)), summary_sha256=digest, field=str(field.relative_to(ROOT)),
                field_sha256=s['artifacts_sha256'][str(field.relative_to(folder))],
                raw_sha256=s['artifacts_sha256'][str(raw_file.relative_to(folder))], velocity_scale=conversion,
                n=n, dx_m=h, liquid_cells=len(ijk), porosity=len(ijk)/n**3,
                carrier_directional_flow_m3_s=target.tolist(),
                native_axial_flow_m3_s=float(target[0]/conversion),
                equivalent_column_flow_m3_s=float(target[0]*inputs['area_m2']/ref.L**2))
    return info, ijk, velocity


def input_arrays(n, h, ijk, velocity):
    # Indices in the handoff are one-based in the original saved fluid row order.
    ids = np.zeros((n,n,n), dtype=np.int32)
    ids[tuple(ijk.T)] = np.arange(1,len(ijk)+1,dtype=np.int32)
    faces = {}
    for axis in range(3):
        destination = np.roll(ids,-1,axis=axis)[tuple(ijk.T)]
        owners = np.flatnonzero(destination)
        neighbors = destination[owners]-1
        flow = .5*h**2*(velocity[owners,axis]+velocity[neighbors,axis])
        faces[axis] = (owners+1, neighbors+1, np.full(len(owners),axis+1), flow)
    arrays = dict(zip(('left.i32','right.i32','axis.u8','unbalanced.f64'),
                     (np.concatenate(column).astype(dtype) for column,dtype in zip(
                         zip(*faces.values(),strict=True), ('<i4','<i4','u1','<f8'),strict=True)),strict=True))
    arrays['varying_initial.f64'] = (.5+ref.A*np.sin(2*np.pi*(ijk[:,0]+.5)/n)).astype('<f8')
    return arrays


def command(directory, info, runtime):
    args = [info['n'],info['dx_m'],ref.D,*info['carrier_directional_flow_m3_s'],*ref.TIMES]
    return [*runtime,str(HERE/'solve.jl'),str(directory),*map(str,args)]


def execute(directory, info, arrays, runtime):
    directory.mkdir(parents=True)
    for name, value in arrays.items():
        value.tofile(directory/name)
    cmd = command(directory, info, runtime)
    write_json(directory/'request.json',dict(command=cmd,source=info,diffusivity_m2_s=ref.D))
    start=perf_counter()
    print('Julia handoff start',info['n'],flush=True)
    with (directory/'stdout.log').open('w') as log:
        p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    write_json(directory/'execution.json',dict(returncode=p.returncode,wall_s=perf_counter()-start))
    if p.returncode:
        raise RuntimeError(f'Julia failed. See {directory / "stdout.log"}')


def flow_metrics(ijk, info, arrays, balanced):
    left, right = arrays['left.i32']-1, arrays['right.i32']-1
    axis, raw = arrays['axis.u8']-1, arrays['unbalanced.f64']
    n, count = info['n'], len(ijk)
    target = np.asarray(info['carrier_directional_flow_m3_s'])
    if balanced.shape != raw.shape or not np.isfinite(balanced).all():
        raise ValueError('Invalid balanced face data')

    def details(flow):
        net = (np.bincount(left,weights=flow,minlength=count)
               -np.bincount(right,weights=flow,minlength=count))
        planes = [np.bincount(ijk[left[axis==d],d],weights=flow[axis==d],minlength=n) for d in range(3)]
        return dict(max_cell_net_flow_m3_s=float(max(abs(net))),
                    rms_cell_net_flow_m3_s=float(np.sqrt(np.mean(net**2))),
                    plane_flows_m3_s=[p.tolist() for p in planes],
                    max_plane_error_relative_to_native_axial=float(np.max(abs(np.array(planes)-target[:,None]))/target[0]))

    before, after = details(raw), details(balanced)
    scale = float(np.max(abs(raw)))
    correction = balanced-raw
    correction_rms = float(np.linalg.norm(correction)/np.linalg.norm(raw))
    return dict(before=before,after=after,face_count=len(raw),max_face_flow_m3_s=scale,
                rms_face_correction_relative=correction_rms,
                max_face_correction_relative_to_max_face=float(max(abs(correction))/scale),
                checks=dict(local_balance=after['max_cell_net_flow_m3_s']<=1e-10*float(max(abs(balanced))),
                    directional_throughput=after['max_plane_error_relative_to_native_axial']<=1e-9,
                    face_fidelity=correction_rms<=.05))


def scalar_metrics(directory, prefix, initial, ijk, n):
    history = np.genfromtxt(directory/(prefix+'_history.csv'),names=True,delimiter=',')
    fields = [np.fromfile(directory/(prefix+'_'+k+'.f64'),dtype='<f8') for k in SNAPSHOTS]
    if any(field.shape != initial.shape or not np.isfinite(field).all() for field in fields):
        raise ValueError('Invalid saved tracer snapshot')
    times = history['time_s']
    if times[0]!=0 or times[-1]!=ref.TIMES[-1] or np.any(np.diff(times)<=0):
        raise ValueError('Invalid tracer time grid')
    if not np.isfinite(history.view(np.float64)).all():
        raise ValueError('Invalid tracer history')
    # Report concentration and measure how much of the initial wave remains.
    drift = float(np.max(abs(history['mean_kg_m3']/initial.mean()-1)))
    lo, hi = float(history['min_kg_m3'].min()), float(history['max_kg_m3'].max())
    amplitude = ref.A if prefix != 'uniform' else 0.
    excursion = max(.5-amplitude-lo,hi-(.5+amplitude),0.)
    checks = dict(mass=drift<=1e-8,engineering_bounds=excursion<=.02,physical_bounds=excursion<=1e-8)
    stats = dict(history_steps=len(times)-1,range_kg_m3=[lo,hi],liquid_inventory_relative_drift=drift,
                 maximum_excursion_kg_m3=excursion,checks=checks)
    phi = 2*np.pi*(ijk[:,0]+.5)/n
    regressors = np.column_stack([np.ones(len(ijk)),np.sin(phi),np.cos(phi)])
    initial_variance = float(np.var(initial))
    stats['snapshots'] = {}
    for key,t,field in zip(SNAPSHOTS,ref.TIMES,fields,strict=True):
        ii = np.flatnonzero(times==t)
        if len(ii)!=1:
            raise ValueError('Missing exact observation time')
        i = ii[0]
        np.testing.assert_allclose([field.mean(),field.min(),field.max()],
            [history['mean_kg_m3'][i],history['min_kg_m3'][i],history['max_kg_m3'][i]],rtol=0,atol=2e-14)
        if amplitude:
            _,a,b = np.linalg.lstsq(regressors,field,rcond=None)[0]
            variance = float(np.var(field))
            diffused = float(history['integrated_dissipation'][i])
            residual = abs(.5*(initial_variance-variance)-diffused)/(.5*initial_variance)
            stats['snapshots'][key] = dict(time_s=t,variance_kg2_m6=variance,variance_retained=variance/initial_variance,
                fundamental_amplitude_kg_m3=float(np.hypot(a,b)),phase_rad=float(np.arctan2(b,a)),
                integrated_dissipation=diffused,variance_budget_relative_residual=residual)
            checks[key+'_variance_budget'] = residual<=1e-5
    if not amplitude:
        checks['constant'] = excursion<=1e-9
    return stats, fields


def evaluate(directory, info, ijk, arrays, runtime):
    req = dict(command=command(directory,info,runtime),source=info,diffusivity_m2_s=ref.D)
    if read_json(directory/'request.json')!=req or read_json(directory/'execution.json')['returncode']:
        raise ValueError('Saved execution request differs')
    for name,value in arrays.items():
        np.testing.assert_array_equal(np.fromfile(directory/name,dtype=value.dtype),value)
    balanced = np.fromfile(directory/'balanced.f64',dtype='<f8')
    n, initial = info['n'], arrays['varying_initial.f64']
    flow = flow_metrics(ijk,info,arrays,balanced)
    uniform,_ = scalar_metrics(directory,'uniform',np.full(len(ijk),.5),ijk,n)
    varying,fields = scalar_metrics(directory,'varying',initial,ijk,n)
    result = dict(source=info,flow=flow,uniform=uniform,varying=varying)
    reference_fields = None
    if n==64:
        reference_fields = [np.fromfile(directory/('exponential_'+k+'.f64'),dtype='<f8') for k in SNAPSHOTS]
        label = 'exponential'
    elif n==96:
        result['half_step'],reference_fields = scalar_metrics(directory,'half_step',initial,ijk,n)
        label = 'time_refinement'
    if reference_fields is not None:
        errors = {k:float(np.sqrt(np.mean(((a-b)/ref.A)**2)))
                  for k,a,b in zip(SNAPSHOTS,fields,reference_fields,strict=True)}
        result[label] = dict(RMS_over_amplitude=errors,passed=all(v<=.001 for v in errors.values()))
    return result


def campaign(base, runtime, *, replay):
    cases = {}
    for n in SOURCES:
        info, ijk, velocity = archived(n)
        directory = base/f'N{n}'
        arrays = input_arrays(n,info['dx_m'],ijk,velocity)
        if not replay:
            execute(directory,info,arrays,runtime)
        cases[str(n)] = evaluate(directory,info,ijk,arrays,runtime)
        if not replay:
            s=cases[str(n)]
            print('Julia handoff complete',n,'correction',s['flow']['rms_face_correction_relative'],
                  'uniform',s['uniform']['range_kg_m3'],'varying',s['varying']['range_kg_m3'],flush=True)
    space = {}
    for time in SNAPSHOTS:
        a,b = (cases[n]['varying']['snapshots'][time] for n in ('64','96'))
        v = abs(a['variance_kg2_m6']/b['variance_kg2_m6']-1)
        m = abs(a['fundamental_amplitude_kg_m3']/b['fundamental_amplitude_kg_m3']-1)
        p = abs(float(np.angle(np.exp(1j*(a['phase_rad']-b['phase_rad'])))))
        space[time] = dict(variance_relative_change=v,amplitude_relative_change=m,phase_change_rad=p,
                           passed=v<=.05 and m<=.05 and p<=.05)
    passed=lambda checks: all(v for k,v in checks.items() if k!='physical_bounds')
    handoff = all(passed(c['flow']['checks']) and passed(c['uniform']['checks']) for c in cases.values())
    spreading = (all(passed(cases[n]['varying']['checks']) for n in ('64','96'))
                 and passed(cases['96']['half_step']['checks']) and cases['64']['exponential']['passed']
                 and cases['96']['time_refinement']['passed'] and all(v['passed'] for v in space.values()))
    return dict(cases=cases,spatial_change=space,handoff_qualified=handoff,spreading_qualified=spreading,
                periodic_control_qualified=handoff and spreading)


def run():
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'runs'/stamp
    base.mkdir(parents=True)
    files=[HERE/k for k in ('run.py','solve.jl','protocol.md')]
    files.extend(ROOT/k for k in ('julia/adsorption/PorousFaceTransport.jl','julia/Project.toml','julia/Manifest.toml'))
    hashes=source_manifest(files)
    for name in hashes:
        target=base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    record=dict(schema='pore2process.porous-face-handoff.v1',run_id=stamp,source_sha256=hashes,
                julia_runtime=ref.runtime(),new_openlb_runs=0,new_bed_runs=0,new_DWSIM_runs=0)
    write_json(base/'plan.json',record)
    print('Porous face-flow handoff',base,flush=True)
    record['analysis']=campaign(base,record['julia_runtime'],replay=False)
    record['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path=base/'summary.json'
    write_json(path,record)
    write_json(OUTPUT/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed porous handoff',path,record['analysis']['periodic_control_qualified'],flush=True)


def check(path):
    path=path.resolve()
    record=validate_record(path)
    result=campaign(path.parent,record['julia_runtime'],replay=True)
    if result!=record['analysis']:
        raise ValueError('Replay differs')
    report=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(record['source_sha256']),
                artifacts=len(record['artifacts_sha256']),handoff_qualified=result['handoff_qualified'],
                spreading_qualified=result['spreading_qualified'])
    write_json(OUTPUT/'verification.json',report)
    print(report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path)
    args=parser.parse_args()
    check(args.check) if args.check else run()
