"""Execute separately guarded stock-native transport qualification cases."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import build_case, verify_vendor
from examples.adsorption.bed import analyze as flow_analysis
from examples.adsorption.packed_bed_tracer.bridge import moments
from .bridge import boundary_moments, reference_profile

HERE=Path(__file__).parent
OUTPUT=ROOT/'results/adsorption/column_transport'


def specifications(hydraulic_run):
    side=.003993871093653167
    jobs=[]
    def add(name,model=0,n=24,D=1e-7,tau=.55,end=2.,initial=1.,periodic=True,
            field='uniform',speed=.0005,scale=1.,outlet='zero',length_side=side):
        jobs.append(dict(id=name,model=model,n=n,D_m2_s=D,tau_diffusive=tau,end_s=end,
            initial_concentration_kg_m3=initial,periodic=periodic,field=str(field),
            speed_m_s=speed,velocity_scale=scale,outlet=outlet,side_m=length_side,cells=1))
    add('constant_free',D=1e-9,tau=.5005)
    for n in (24,32,48):
        add(f'constant_fcc_grid{n}',n=n,field=ROOT/f'results/adsorption/bed/grid{n}/field.csv')
    add('constant_fcc_time24',tau=.525,field=ROOT/'results/adsorption/bed/grid24/field.csv')
    for model,label in [(0,'trt27'),(1,'regularized27'),(2,'mrt7'),(3,'trt7')]:
        add('molecular_fcc_'+label,model=model,D=1e-9,tau=.5005,
            field=ROOT/'results/adsorption/bed/grid24/field.csv')
    for model,outlet in [(0,'zero'),(3,'robin')]:
        add('outlet_constant_'+outlet,model=model,D=1e-7,tau=.55,end=40.,periodic=False,outlet=outlet)
        for n in (16,32):
            add(f'outlet_step_{outlet}{n}',model=model,n=n,D=1e-7,tau=.55,end=40.,initial=0.,periodic=False,outlet=outlet)
    # The hydraulic source remains an actual separately converged field. Its
    # carrier scaling preserves Re; no extrapolated fluid profile is accepted.
    p=Path(hydraulic_run)
    carrier=json.loads((p/'inputs.json').read_text())
    candidates=sorted((p/'native_flow/d800_grid32').glob('iteration*/raw.json'),key=lambda x:int(x.parent.name.removeprefix('iteration')))
    source=candidates[-1].parent
    flow=flow_analysis(source)
    if abs(flow['superficial_velocity_m_s']/carrier['native_target_velocity_m_s']-1)>1e-4:
        raise ValueError('Need the completed fixed-flow d800 N32 native source.')
    physical_scale=carrier['carrier_viscosity_Pa_s']/carrier['carrier_density_kg_m3']/1e-6
    velocity=np.genfromtxt(source/'field.csv',delimiter=',',names=True)
    maxu=float(np.max(np.sqrt(sum(velocity[k]**2 for k in ('ux_m_s','uy_m_s','uz_m_s')))))*physical_scale
    side=flow['raw']['unit_cell_side_m'];h=side/32
    dt=.02*h/maxu
    v=flow['superficial_velocity_m_s']*physical_scale/flow['voxel_porosity']
    for model,label in [(0,'trt27'),(1,'regularized27'),(2,'mrt7'),(3,'trt7')]:
        invcs=4 if model in (2,3) else 3
        add('application_constant_'+label,model=model,n=32,D=1e-9,
            tau=.5+1e-9*dt*invcs/h**2,end=2*side/v,
            field=source/'field.csv',speed=v,scale=physical_scale,length_side=side)
    return jobs,dict(hydraulic_run=str(p.relative_to(ROOT)),field_sha256=sha(source/'field.csv'),
        flow_raw_sha256=sha(source/'raw.json'),physical_velocity_scale=physical_scale,
        application_particle_Peclet=v*.0008/1e-9,
        application_max_cell_Peclet=maxu*h/1e-9,
        source_role='Synthetic FCC at source flow/diameter; Dm=1e-9 is a declared phenol-water scenario.')


def execute(base,spec,binary):
    dest=base/'native'/spec['id'];dest.mkdir(parents=True)
    command=[str(binary),str(spec['n']),str(spec['cells']),str(spec['side_m']),str(spec['D_m2_s']),
        str(spec['tau_diffusive']),str(spec['end_s']),spec['field'],str(spec['speed_m_s']),
        str(spec['initial_concentration_kg_m3']),str(spec['velocity_scale']),str(int(spec['periodic'])),spec['outlet']]
    write_json(dest/'request.json',dict(spec=spec,command=command,binary_sha256=sha(binary),
        native_field_sha256=None if spec['field']=='uniform' else sha(Path(spec['field']))))
    print('Native start',spec['id'],flush=True);start=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=1800)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-start))
    if result.returncode:raise RuntimeError('Native execution failed: '+str(dest/'stdout.log'))
    print('Native finish',spec['id'],f'{perf_counter()-start:.1f}s',flush=True)


def analyze(path, *, references=True):
    path=Path(path);raw=json.loads((path/'raw.json').read_text())
    a=np.genfromtxt(path/'inventory.csv',delimiter=',',names=True)
    f=np.genfromtxt(path/'flux.csv',delimiter=',',names=True)
    scale=max(raw['initial_inventory_kg'],float(a['cumulative_inlet_kg'][-1]),1e-30)
    residual=a['fluid_tracer_kg']+a['wall_population_tracer_kg']-raw['initial_inventory_kg']-a['cumulative_inlet_kg']+a['cumulative_outlet_kg']
    np.testing.assert_allclose(residual,a['mass_residual_kg'],atol=1e-20,rtol=1e-8)
    np.testing.assert_allclose(np.cumsum(f['inlet_flux_kg_s'])*raw['dt_s'],a['cumulative_inlet_kg'][1:],atol=1e-18,rtol=1e-9)
    np.testing.assert_allclose(np.cumsum(f['outlet_flux_kg_s'])*raw['dt_s'],a['cumulative_outlet_kg'][1:],atol=1e-18,rtol=1e-9)
    finite=bool(all(np.all(np.isfinite(a[k])) for k in a.dtype.names))
    low=float(np.min(a['min_concentration_kg_m3']));high=float(np.max(a['max_concentration_kg_m3']))
    constant=raw['initial_concentration']==1.
    departure=max(abs(low-1),abs(high-1)) if constant else None
    checks=dict(finite=finite,mass=float(np.max(abs(residual)))/scale<=1e-8,
        positive=low>=-1e-8,bounded=high<=1+1e-6,
        uniform_preservation=None if not constant else departure<=1e-6)
    r=dict(raw=raw,checks=checks,minimum_concentration_kg_m3=low,maximum_concentration_kg_m3=high,
        uniform_maximum_departure=departure,normalized_mass_residual=float(np.max(abs(residual)))/scale,
        role='Constant preservation counterexample' if constant else 'Stock outlet step reference',
        qualification_passed=all(v is not False for v in checks.values()))
    if not constant and raw['uniform'] and references:
        m=moments(f['time_s'].tolist(),f['outlet_flux_kg_s'].tolist())
        L=(raw['nx']-2)*raw['dx_m'];v=raw['flow_m3_s']/raw['area_m2']
        p=dict(length_m=L,velocity_m_s=v,dispersion_m2_s=raw['diffusivity_m2_s'])
        exact=boundary_moments('dirichlet_neumann',**p)
        model=reference_profile(**p,t_end_s=float(a['time_s'][-1]),dt_s=.002,cells=800)
        pred=np.interp(f['time_s'],model['time_s'],np.array(model['outlet_flux_per_area_m_s'])/v)
        y=f['outlet_flux_kg_s']/raw['flow_m3_s']
        comparisons=dict(mean_relative_error=abs(m['mean_s']/exact['mean_s']-1),
            variance_relative_error=abs(m['variance_s2']/exact['variance_s2']-1),
            profile_RMSE=float(np.sqrt(np.mean((pred-y)**2))),
            absolute_plateau_relative_error=abs(m['steady_exit_flux_kg_s']/raw['flow_m3_s']-1))
        checks.update(reference_mean=comparisons['mean_relative_error']<=.05,
            reference_variance=comparisons['variance_relative_error']<=.1,
            reference_profile=comparisons['profile_RMSE']<=.02,
            absolute_plateau=comparisons['absolute_plateau_relative_error']<=.01,
            tail=m['tail_second_moment_fraction']<.01 and m['tail_area_fraction']<.01,
            final_signal=m['final_fraction']>=.999)
        r.update(moments=m,analytical=exact,comparison=comparisons,reference_inputs=p,
            reference_normalized_mass_residual=model['normalized_mass_residual'])
        r['qualification_passed']=all(v is not False for v in checks.values())
    return r


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--hydraulic-run',type=Path,required=True)
    args=parser.parse_args()
    start=perf_counter();rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'transport/runs'/rid;base.mkdir(parents=True)
    jobs,application=specifications(args.hydraulic_run.resolve())
    files=[Path(__file__),HERE/'bridge.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir()),
        ROOT/'julia/adsorption/ColumnBoundaryReference.jl',ROOT/'julia/adsorption/test_column_boundary_reference.jl',
        ROOT/'examples/adsorption/packed_bed_tracer/bridge.py',ROOT/'julia/adsorption/TracerContact.jl']
    sources=source_manifest(files)
    for name in sources:
        target=base/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    binaries={};builds={}
    for model in range(4):
        src=base/f'generated_source/model{model}';shutil.copytree(HERE/'openlb',src)
        config=src/'config.mk';config.write_text(config.read_text()+f'\nCXXFLAGS += -DP2P_SCALAR_MODEL={model}\n')
        build=ROOT/f'.local/adsorption/column_transport_model{model}'
        builds[str(model)]=build_case(src,build,base/f'build{model}')
        binaries[model]=build/'capillary3d'
    write_json(base/'application.json',application)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda s:execute(base,s,binaries[s['model']]),jobs))
    cases={s['id']:analyze(base/'native'/s['id']) for s in jobs}
    packed=[c for key,c in cases.items() if 'fcc' in key or key.startswith('application_')]
    qualified=all(c['qualification_passed'] for c in packed)
    closure=dict(boundary='convective-column',native_run_id=rid,
        calibration_boundary='native-inert-tracer-only',axial_dispersion_m2_s=None,
        checks=dict(uniform_preservation=qualified,mass=all(c['checks']['mass'] for c in cases.values()),
            bounds=qualified,mesh=False,time=False,convective_outlet=False,longer_segment_prediction=False),
        reason='No transport coefficient is extracted when uniform preservation fails.')
    s=dict(schema='pore2process.column-transport-audit.v1',run_id=rid,source_sha256=sources,
        cases=cases,application=application,closure=closure,transport_qualified=qualified,
        native_runs=len(jobs),new_DWSIM_runs=0,new_adsorption_predictions=0,
        builds=builds,vendor=verify_vendor(),elapsed_s=perf_counter()-start)
    s['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',s)
    write_json(OUTPUT/'transport/latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Transport saved',base,'qualified',qualified,flush=True)
    print({k:v['checks'] for k,v in cases.items()},flush=True)


if __name__=='__main__':main()
