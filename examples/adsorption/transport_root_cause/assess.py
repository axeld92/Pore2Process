"""Replay guarded mesh evidence and compare exact discrete flux diagnostics."""
import argparse
import json
from pathlib import Path

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from .study import HERE, OUTPUT, check
from .flux_audit.analysis import summarize, summarize_uniform_stencil
from .flux_audit.findings import peak_statistics


def peak_departure(path):
    """Measure affected liquid volume and the magnitude beyond the bounds."""
    f=np.genfromtxt(path,delimiter=',',names=True,usecols=('material','concentration_kg_m3'))
    c=f['concentration_kg_m3'][f['material']!=2]
    if not c.size or not np.all(np.isfinite(c)) or np.sum(c)<=0:
        raise ValueError('Peak needs finite fluid concentration and positive liquid inventory')
    above=np.maximum(c-.75,0.);below=np.maximum(.25-c,0.)
    inventory=float(np.sum(c))
    # Equal lattice cell volumes cancel in these dimensionless ratios.
    return dict(fluid_nodes=len(c),fraction_liquid_beyond_engineering_allowance=float(np.mean((above>.02)|(below>.02))),
        above_physical_maximum_over_liquid_inventory=float(np.sum(above))/inventory,
        below_physical_minimum_over_liquid_inventory=float(np.sum(below))/inventory,
        absolute_beyond_bounds_over_liquid_inventory=float(np.sum(above+below))/inventory,
        note='Selected peak field only; spatial error indicator, not mass loss or a process-accuracy bound.')


def history_departure(path):
    h=np.genfromtxt(path,delimiter=',',names=True)
    error=np.maximum(0.,np.maximum(.25-h['min_concentration_kg_m3'],h['max_concentration_kg_m3']-.75))
    violating=np.flatnonzero(error>.02)
    return dict(first_engineering_violation_time_s=None if not len(violating) else float(h['time_s'][violating[0]]),
        fraction_sampled_steps_outside_engineering_screen=float(len(violating)/len(h)),
        final_excursion_kg_m3=float(error[-1]),maximum_excursion_kg_m3=float(np.max(error)))


def assess(paths):
    records=[];cases={}
    for path in paths:
        path=Path(path).resolve();s=check(path)
        records.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),
                            native_runs=s['new_native_runs'],sources=len(s['source_sha256']),artifacts=len(s['artifacts_sha256'])))
        for name, result in s['cases'].items():
            if name in cases:raise ValueError('Repeated case ID: '+name)
            dest=path.parent/'native'/name
            r=result['raw'];audit_files=sorted(dest.glob('audit_*.csv'))
            audits={f.stem:summarize(f,dx_m=r['dx_m'],dt_s=r['dt_s']) for f in audit_files}
            constant=result['peak']['initial_amplitude_kg_m3']==0
            stencil=summarize_uniform_stencil(dest/'uniform_stencil.csv',
                first_step_audit=dest/'audit_1.csv' if constant else None,initial_concentration=.5)
            cases[name]=dict(result=result,exact_population_audits=audits,uniform_stencil=stencil,
                             native_directory=str(dest.relative_to(ROOT)))
    mesh={k:v for k,v in cases.items() if k in ('n32_sine','n48_sine','n64_sine','n96_sine')}
    for item in mesh.values():
        base=ROOT/item['native_directory']
        item['peak_departure']=peak_departure(base/'peak.csv')
        item['history_departure']=history_departure(base/'inventory.csv')
        if abs(item['history_departure']['maximum_excursion_kg_m3']-item['result']['peak']['maximum_excursion_kg_m3'])>1e-12:
            raise ValueError('Saved peak disagrees with all-time history')
    if 'n96_sine' in mesh:
        fine=mesh['n96_sine'];raw=fine['result']['raw']
        fine['peak_geometry']=peak_statistics(ROOT/fine['native_directory']/'peak.csv',
            side_m=raw['side_m'],particle_radius_m=.0004)
    comparisons=[]
    keys=sorted(mesh,key=lambda name:int(name.split('_')[0][1:]))
    for a,b in zip(keys,keys[1:]):
        x,y=mesh[a]['result'],mesh[b]['result']
        ex,ey=x['peak']['maximum_excursion_kg_m3'],y['peak']['maximum_excursion_kg_m3']
        comparisons.append(dict(coarse=a,fine=b,excursion_reduction_fraction=None if ex==0 else 1-ey/ex,
            time_step_relative_difference=abs(x['raw']['dt_s']/y['raw']['dt_s']-1),
            voxel_porosity_coarse=x['raw']['fluid_volume_m3']/x['raw']['side_m']**3,
            voxel_porosity_fine=y['raw']['fluid_volume_m3']/y['raw']['side_m']**3,
            liquid_SD_relative_change=y['periodic_wave_final']['sd_C_kg_m3']/x['periodic_wave_final']['sd_C_kg_m3']-1,
            first_harmonic_projection_norm_relative_change=np.hypot(y['periodic_wave_final']['sin_coefficient_kg_m3'],y['periodic_wave_final']['cos_coefficient_kg_m3'])/np.hypot(x['periodic_wave_final']['sin_coefficient_kg_m3'],x['periodic_wave_final']['cos_coefficient_kg_m3'])-1))
    # Compare added instrumentation with preserved, otherwise identical runs.
    parent=json.loads((ROOT/records[0]['path']).read_text())
    old=json.loads((ROOT/parent['parent_summary']).read_text())
    repeat={}
    for name, old_name in [('n32_sine','half_dt'),('n48_sine','grid48_matched_dt')]:
        if name not in cases:continue
        a,b=cases[name]['result'],old['cases'][old_name]
        repeat[name]={k:abs(a[k]-b[k]) for k in ('minimum_concentration_kg_m3','maximum_concentration_kg_m3')}
    result=dict(schema='pore2process.transport-root-cause-assessment.v1',parents=records,cases=cases,
                analysis_source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                    (Path(__file__),HERE/'study.py',HERE/'flux_audit/analysis.py',HERE/'flux_audit/findings.py')},
                mesh_comparisons=comparisons,instrumentation_reproducibility=repeat,
                actual_new_native_runs=sum(x['native_runs'] for x in records),new_native_runs_in_check=0,
                new_DWSIM_runs=0,new_adsorption_predictions=0,service_prediction_qualified=False,
                qualification_note='Short periodic diagnostics cannot qualify failed open-boundary transport or supply column dispersion.')
    write_json(OUTPUT/'assessment.json',result)
    write_json(HERE/'reference/investigation_reference.json',result)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('summaries',nargs='*',type=Path)
    p.add_argument('--manifest',action='store_true');args=p.parse_args()
    paths=args.summaries
    if args.manifest:
        if paths:p.error('Choose explicit summaries or --manifest')
        manifest=json.loads((HERE/'reference/study_manifest.json').read_text())
        paths=[]
        for entry in manifest['scalar_summaries']:
            path=ROOT/entry['path']
            if sha(path)!=entry['sha256']:raise ValueError('Manifest scalar summary changed')
            paths.append(path)
    if not paths:p.error('Provide summaries or --manifest')
    r=assess(paths)
    print(json.dumps(dict(mesh_comparisons=r['mesh_comparisons'],instrumentation_reproducibility=r['instrumentation_reproducibility'],
        actual_new_native_runs=r['actual_new_native_runs'],new_native_runs_in_check=0),indent=2))


if __name__=='__main__':main()
