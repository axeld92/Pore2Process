"""Verify the complete diagnostic dependency graph without new native runs."""
import json
from pathlib import Path

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha, write_json
from pore2process_demo.provenance import verified_source_path
from pore2process_demo.runtime.openlb import verify_vendor
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.bed import analyze as analyze_flow
from examples.adsorption.transport_accuracy.diagnosis import field_audit
from .study import HERE, OUTPUT, check as check_scalar
from .hydraulics.run import native_checks
from .flux_audit.analysis import summarize


def verify_dependencies(path, visited):
    path=Path(path).resolve()
    if path in visited:return visited[path]
    s=validate_record(path);visited[path]=s
    for role in ('parent','hydraulic','refined','serial','native_crosscheck','previous_refinement','baseline','flow'):
        key=role+'_summary'
        if not s.get(key):continue
        digest_key='native_crosscheck_sha256' if role=='native_crosscheck' else key+'_sha256'
        if role=='previous_refinement':digest_key='previous_refinement_sha256'
        p=ROOT/s[key]
        if sha(p)!=s[digest_key]:raise ValueError('Changed dependency: '+str(p))
        verify_dependencies(p,visited)
    for entry in s.get('extra_summaries',[]):
        p=ROOT/entry['path']
        if sha(p)!=entry['sha256']:raise ValueError('Changed explicit extra dependency')
        verify_dependencies(p,visited)
    # The flat-wall continuation preserves a failed-oracle cohort through an
    # explicit manifest instead of pretending that cohort had a passed summary.
    if s.get('parent_manifest'):
        p=path.parent/s['parent_manifest']
        if sha(p)!=s['parent_manifest_sha256']:raise ValueError('Changed parent artifact manifest')
        parent=json.loads(p.read_text());base=ROOT/parent['parent']
        for name,digest in parent['source_sha256'].items():
            if sha(verified_source_path(name,digest))!=digest:raise ValueError('Changed manifest parent source')
        for name,digest in parent['artifacts_sha256'].items():
            if sha(base/name)!=digest:raise ValueError('Changed manifest parent artifact')
        visited[p.resolve()]=parent
    return s


def verify_scalar(path, visited):
    path=Path(path).resolve();s=check_scalar(path);verify_dependencies(path,visited)
    if s.get('schema')=='pore2process.matched-constant-background.v1':
        from .constant import comparison
        native=path.parent/'native'/s['specifications'][0]['id']
        if comparison(ROOT/s['parent_summary'],native)!=s['background_comparison']:
            raise ValueError('Matched-background diagnostic does not replay')
    for spec in s['specifications']:
        audit=field_audit(spec,path.parent/'native'/spec['id']/'final.csv')
        if max(audit['coordinate_maximum_error_m'].values())>1e-15:
            raise ValueError('Native/scalar physical coordinates differ')
        r=s['cases'][spec['id']]['raw']
        if r['max_lattice_speed']>.05:
            raise ValueError('Scalar lattice speed exceeds the declared diagnostic range')
        recovered_D=(r['magic']/(r['tau_even']-.5))*r['dx_m']**2/(3*r['dt_s'])
        if abs(recovered_D/1e-9-1)>1e-9:
            raise ValueError('Scalar physical diffusivity changed')
    if s.get('source_point_record'):
        if sha(ROOT/s['source_point_record'])!=s['source_point_record_sha256']:
            raise ValueError('Changed source-point record')
    if s.get('refined_summary'):
        refined=visited[(ROOT/s['refined_summary']).resolve()]
        point=json.loads((path.parent/'source_point.json').read_text())
        if point!=refined['source_point']:raise ValueError('Scalar point differs from accepted hydraulic summary')
        carrier=visited[(ROOT/s['hydraulic_summary']).resolve()]['inputs']
        for spec in s['specifications']:
            for key,value in point.items():
                if spec[key]!=value:raise ValueError('Scalar specification changed source-point '+key)
            if sha(Path(spec['field']))!=point['field_sha256']:raise ValueError('Changed hydraulic field')
            if sha(Path(spec['field']).parent/'raw.json')!=point['raw_sha256']:raise ValueError('Changed hydraulic raw result')
            if not all(native_checks(analyze_flow(Path(spec['field']).parent),carrier['native_target_velocity_m_s']).values()):
                raise ValueError('Unqualified fixed-flow source')
    if s.get('serial_summary'):
        serial_path=ROOT/s['serial_summary']
        for spec in s['specifications']:
            old=serial_path.parent/'native'/f"n{spec['n']}_sine"
            new=path.parent/'native'/spec['id'];comparisons={}
            for name in ('peak.csv','final.csv','inventory.csv','moments.csv'):
                a=np.genfromtxt(old/name,delimiter=',',names=True);b=np.genfromtxt(new/name,delimiter=',',names=True)
                if a.shape!=b.shape or a.dtype.names!=b.dtype.names:raise ValueError('Parallel equivalence shape changed')
                for field in a.dtype.names:np.testing.assert_allclose(a[field],b[field],rtol=1e-12,atol=2e-12)
                comparisons[name]=max(float(np.max(abs(a[field]-b[field]))) for field in a.dtype.names)
            if comparisons!=s['serial_equivalence_maximum_absolute_differences']:raise ValueError('Equivalence does not replay')
    return s


def main():
    manifest_path=HERE/'reference/study_manifest.json';manifest=json.loads(manifest_path.read_text())
    visited={};scalar=[]
    for entry in manifest['scalar_summaries']:
        p=ROOT/entry['path']
        if sha(p)!=entry['sha256']:raise ValueError('Manifest scalar changed')
        scalar.append(verify_scalar(p,visited))
    water={}
    for entry in manifest.get('hydraulic_summaries',[]):
        p=ROOT/entry['path']
        if sha(p)!=entry['sha256']:raise ValueError('Manifest hydraulic changed')
        s=verify_dependencies(p,visited)
        if not s['checks_passed']:raise ValueError('Hydraulic source failed its screens')
        # Audit every completed pressure-root iteration, with its actual dt.
        for iteration in s['native']['iterations']:
            dest=p.parent/iteration['directory'];r=iteration['analysis']['raw']
            f=dest/'water_population_audit.csv'
            if not f.exists():
                f=next(dest.glob('*population*.csv'))
            a=summarize(f,dx_m=r['dx_m'],dt_s=r['dt_s'],density_scale=1000.)
            if not all(a['checks'].values()):raise ValueError('Water population accounting failed')
            water[str(dest.relative_to(ROOT))]=a
    for role in ('channel_summary','thin_summary','candidate_summary','crosscheck_summary','interrupted_summary','flat_summary','phase_summary','halfway_summary','halfway_review'):
        entry=manifest.get(role)
        if not entry:continue
        p=ROOT/entry['path']
        if sha(p)!=entry['sha256']:raise ValueError('Changed '+role)
        verify_dependencies(p,visited)
    for entry in manifest.get('compatibility_summaries',[]):
        p=ROOT/entry['path']
        if sha(p)!=entry['sha256']:raise ValueError('Changed compatibility assessment')
        verify_dependencies(p,visited)
    for path,s in list(visited.items()):
        if s.get('status')=='interrupted_by_host_session_update':
            for iteration in s['native']['iterations']:
                dest=path.parent/iteration['directory'];r=iteration['analysis']['raw']
                a=summarize(dest/'water_population_audit.csv',dx_m=r['dx_m'],dt_s=r['dt_s'],density_scale=1000.)
                if not all(a['checks'].values()):raise ValueError('Preserved cold water population audit failed')
                water[str(dest.relative_to(ROOT))]=a
        if s.get('schema')!='pore2process.stock-openmp-crosscheck.v1':continue
        for key,case in s['cases'].items():
            dest=path.parent/key;r=case['analysis']['raw']
            f=dest/'water_population_audit.csv'
            a=summarize(f,dx_m=r['dx_m'],dt_s=r['dt_s'],density_scale=1000.)
            if not all(a['checks'].values()):raise ValueError('Water cross-check population audit failed')
            water[str(dest.relative_to(ROOT))]=a
    sources={};artifacts=0
    for path,s in visited.items():
        sources.update(s['source_sha256']);artifacts+=len(s['artifacts_sha256'])
    cases={k:v for cohort in scalar for k,v in cohort['cases'].items()}
    mesh_names=[k for k in ('n32_sine','n48_sine','n64_sine','n96_sine') if k in cases]
    raw=[cases[k]['raw'] for k in mesh_names]
    dt_error=max(abs(r['dt_s']/raw[0]['dt_s']-1) for r in raw)
    duration_error=max(abs(r['dt_s']*r['steps']/(raw[0]['dt_s']*raw[0]['steps'])-1) for r in raw)
    if max(dt_error,duration_error)>1e-9:raise ValueError('Matched spatial series changed time step or window')
    mesh_check=dict(cases=mesh_names,maximum_relative_time_step_difference=dt_error,
        maximum_relative_duration_difference=duration_error)
    if 'n96_sine' in cases:
        binaries={}
        for entry,cohort in zip(manifest['scalar_summaries'],scalar):
            for name in set(cohort['cases'])&{'n64_sine','n96_sine'}:
                request=ROOT/Path(entry['path']).parent/'native'/name/'request.json'
                binaries[name]=json.loads(request.read_text())['binary_sha256']
        if binaries['n96_sine']!=binaries['n64_sine']:raise ValueError('Fine mesh changed scalar executable')
        mesh_check['N64_and_N96_executable_sha256']=binaries['n96_sine']
    policy=manifest.get('conditional_fine_time_check');time_check=None
    if policy:
        if sha(ROOT/policy['path'])!=policy['sha256']:raise ValueError('Changed prospective fine time-step decision')
        if 'n96_sine' in cases:
            excursion=cases['n96_sine']['peak']['maximum_excursion_kg_m3']
            trigger=max(policy['engineering_excursion_limit_kg_m3'],policy['strong_improvement_trigger_kg_m3'])
            required=excursion<=trigger;performed='n96_sine_half_dt' in cases
            time_check=dict(observed_excursion_kg_m3=excursion,trigger_kg_m3=trigger,required=required,performed=performed)
            if manifest['status']=='complete' and required!=performed:
                raise ValueError('Fine time-step result differs from the prospective decision')
        elif manifest['status']=='complete':raise ValueError('Completed manifest requires a matched N96 scalar run')
    result=dict(schema='pore2process.transport-root-cause-verification.v1',manifest_sha256=sha(manifest_path),
        checker_sha256=sha(Path(__file__)),
        guarded_records=len(visited),unique_source_files_checked=len(sources),artifact_checks=artifacts,
        scalar_native_runs=sum(s['new_native_runs'] for s in scalar),water_population_audits=water,
        matched_spatial_series=mesh_check,fine_time_check=time_check,
        new_native_runs_in_check=0,new_DWSIM_runs_in_check=0,vendor=verify_vendor(),
        service_prediction_qualified=False)
    write_json(OUTPUT/'verification.json',result)
    print('Verified',len(visited),'guarded records,',len(sources),'sources,',artifacts,'artifact checks; 0 new native/DWSIM runs',flush=True)


if __name__=='__main__':main()
