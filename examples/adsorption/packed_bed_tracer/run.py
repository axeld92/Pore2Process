"""Execute native cases separately from Julia translation and saved-data checks."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest
from pore2process_demo.io import write_json,sha
from pore2process_demo.runtime.openlb import build_case,verify_vendor
from pore2process_demo.bridge import runtime_versions
from .analyze import analyze
from .bridge import call
from examples.adsorption.bed import analyze as analyze_flow

HERE=Path(__file__).parent
BUILD=ROOT/'.local/adsorption/packed_bed_tracer'


def specs(flow16):
    bed=json.loads((ROOT/'results/adsorption/bed/summary.json').read_text())
    if not all(bed['checks'].values()):raise ValueError('Prerequisite native bed checks failed.')
    jobs=[]
    bed['cases']['grid16']=analyze_flow(flow16)
    for n in (16,32):
        r=bed['cases'][f'grid{n}']['raw'];side=r['unit_cell_side_m']
        jobs.append(dict(id=f'uniform{n}',n=n,cells=1,side_m=side,Dm=2e-6,tau=1.,end_s=35.,
            field='uniform',speed=.0005,flow_source=None))
    for name,n,cells,tau,Dm in [('bed16',16,1,1.,2e-6),('bed24',24,1,1.,2e-6),
        ('bed32',32,1,1.,2e-6),('time32',32,1,.75,2e-6),('length32',32,2,1.,2e-6),('diffusion24',24,1,1.,4e-6)]:
        src=flow16 if n==16 else ROOT/f'results/adsorption/bed/grid{n}'
        r=bed['cases'][f'grid{n}'];side=r['raw']['unit_cell_side_m']
        assert r['raw']['converged'] and r['raw']['collision']=='TRT'
        jobs.append(dict(id=name,n=n,cells=cells,side_m=side,Dm=Dm,tau=tau,
            end_s=35.*cells,
            field=str(src/'field.csv'),speed=r['superficial_velocity_m_s']/r['voxel_porosity'],
            flow_source=dict(raw_path=str((src/'raw.json').relative_to(ROOT)),raw_sha256=sha(src/'raw.json'),
                field_sha256=sha(src/'field.csv'),source_summary_sha256=sha(ROOT/'results/adsorption/bed/summary.json'))))
    return jobs


def execute(base,spec):
    dest=base/'native'/spec['id'];dest.mkdir(parents=True)
    command=[str(BUILD/'capillary3d'),str(spec['n']),str(spec['cells']),str(spec['side_m']),str(spec['Dm']),
        str(spec['tau']),str(spec['end_s']),spec['field'],str(spec['speed'])]
    write_json(dest/'request.json',dict(spec=spec,command=command,binary_sha256=sha(BUILD/'capillary3d')))
    print('Starting native',spec['id'],flush=True)
    started=perf_counter()
    with (dest/'stdout.log').open('w') as stream:
        result=subprocess.run(command,cwd=dest,stdout=stream,stderr=subprocess.STDOUT,timeout=2400)
    write_json(dest/'execution.json',dict(returncode=result.returncode,elapsed_s=perf_counter()-started))
    if result.returncode:raise RuntimeError(f'Native {spec["id"]} failed: {dest/"stdout.log"}')
    print('Finished native',spec['id'],f'{perf_counter()-started:.1f} s',flush=True)
    return spec['id']


def relative(a,b):
    return {key:abs(a['moments'][key]/b['moments'][key]-1) for key in ('mean_s','variance_s2')}


def translate(base,cases):
    # Only the short highest-resolution native segment determines coefficients.
    a=cases['bed32'];r=a['raw'];m=a['moments']
    contact=call('derive_contact',length_m=a['transport_length_m'],mean_s=m['mean_s'],variance_s2=m['variance_s2'])
    contact['calibration_case']='bed32';contact['holdout_case']='length32'
    contact['role']='effective contact coefficients from synthetic inert tracer with specified absorbing boundaries; not experimental fluid dispersion'
    predictions={}
    for key in ('bed32','length32'):
        raw=cases[key]['raw'];flux=np.genfromtxt(base/'native'/key/'flux.csv',delimiter=',',names=True)
        p=dict(length_m=cases[key]['transport_length_m'],velocity_m_s=contact['effective_velocity_m_s'],
            dispersion_m2_s=contact['effective_dispersion_m2_s'],t_end_s=float(flux['time_s'][-1]),dt_s=.002,cells=800)
        model=call('solve_ade',**p)
        exact=call('absorbing_moments',length_m=p['length_m'],velocity_m_s=p['velocity_m_s'],dispersion_m2_s=p['dispersion_m2_s'])
        response=np.interp(flux['time_s'],model['time_s'],np.array(model['outlet_flux_per_area_m_s'])/exact['steady_flux_per_area_m_s'])
        native=flux['outlet_flux_kg_s']/cases[key]['moments']['steady_exit_flux_kg_s']
        metrics=dict(RMSE=float(np.sqrt(np.mean((response-native)**2))),
            mean_relative_error=abs(exact['mean_s']/cases[key]['moments']['mean_s']-1),
            variance_relative_error=abs(exact['variance_s2']/cases[key]['moments']['variance_s2']-1))
        metrics['passed']=metrics['RMSE']<=.03 and metrics['mean_relative_error']<=.1 and metrics['variance_relative_error']<=.2
        predictions[key]=dict(inputs=p,result=model,metrics=metrics,analytical_moments=exact)
        write_json(base/f'julia_{key}.json',predictions[key])
    contact['physical_void_velocity_m_s']=r['flow_m3_s']/r['fluid_volume_m3']*((r['nx']-2)*r['dx_m'])
    contact['effective_to_physical_void_velocity_ratio']=contact['effective_velocity_m_s']/contact['physical_void_velocity_m_s']
    contact['effective_contact_volume_m3']=r['flow_m3_s']*m['mean_s']
    contact['physical_fluid_volume_m3']=r['fluid_volume_m3']
    contact['inventory_substitution_applied']=False
    return contact,predictions


def main():
    started=perf_counter();rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=ROOT/'results/adsorption/packed_bed_tracer/runs'/rid;base.mkdir(parents=True)
    build=build_case(HERE/'openlb',BUILD,base/'build')
    flow_build=ROOT/'.local/adsorption/packed_bed_tracer_flow'
    flow_metadata=build_case(ROOT/'examples/adsorption/openlb_bed',flow_build,base/'flow_build')
    flow16=base/'flow16';flow16.mkdir()
    command=[str(flow_build/'capillary3d'),'16','10','.5144','.0023']
    write_json(flow16/'request.json',dict(command=command,binary_sha256=sha(flow_build/'capillary3d')))
    print('Starting native flow16',flush=True)
    with (flow16/'stdout.log').open('w') as stream:
        subprocess.run(command,cwd=flow16,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=2400)
    f=analyze_flow(flow16)
    if f['plane_mass_relative_spread']>1e-6 or max(f['stationarity'].values())>1.0001e-7:
        raise ValueError('New coarse native flow conservation/stationarity failed.')
    jobs=specs(flow16)
    files=[*HERE.glob('*.py'),HERE/'protocol.md',*list((HERE/'openlb').iterdir()),
        ROOT/'julia/adsorption/TracerContact.jl',ROOT/'julia/adsorption/test_tracer_contact.jl',
        ROOT/'tests/adsorption/test_packed_bed_tracer.py',ROOT/'examples/adsorption/bed.py',
        *list((ROOT/'examples/adsorption/openlb_bed').iterdir()),
        ROOT/'uv.lock',ROOT/'julia/Project.toml',ROOT/'julia/Manifest.toml']
    sources=source_manifest(files)
    for name in sources:
        target=base/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,target)
    write_json(base/'source_sha256.json',sources)
    # Native subprocesses may run independently. Julia is invoked only on this
    # main thread AFTER all native jobs finish; no JuliaCall background tasks.
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda job:execute(base,job),jobs))
    cold=perf_counter();cases={spec['id']:analyze(base/'native'/spec['id']) for spec in jobs};cold=perf_counter()-cold
    refinements=dict(grid=relative(cases['bed32'],cases['bed24']),time=relative(cases['time32'],cases['bed32']))
    native_ok=all(all(c['checks'].values()) for c in cases.values())
    numerical=all(v['mean_s']<=.05 and v['variance_s2']<=.1 for v in refinements.values())
    contact,predictions=translate(base,cases)
    summary=dict(schema='pore2process.packed-bed-tracer.v1',run_id=rid,source_sha256=sources,cases=cases,
        refinements=refinements,native_checks_passed=native_ok,numerical_refinement_passed=numerical,
        contact=contact,macro_comparison={key:r['metrics'] for key,r in predictions.items()},
        process_scope='No new DWSIM solve or adsorption calibration; native inert-tracer contact boundary only.',
        new_native_scalar_runs=len(jobs),new_native_flow_runs=1,new_DWSIM_runs=0,
        build=build,flow_build=flow_metadata,new_flow16=f,vendor=verify_vendor(),versions=runtime_versions(),analysis_and_initialization_s=cold,
        total_s=perf_counter()-started)
    summary['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',summary)
    write_json(base.parent.parent/'latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('Saved',base,'native',native_ok,'refinement',numerical,'macro',summary['macro_comparison'],flush=True)

if __name__=='__main__':main()
