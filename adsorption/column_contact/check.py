"""Replay the actual Julia model and reconstruct inventories and selection."""
import json
import shutil
import numpy as np
from adsorption.alhamed.source import ROOT,sha
from adsorption.alhamed.study import observations,compare
from adsorption.film_particle.check import verify_hashes,inventories,close
from .audit import audit,frozen
from .bridge import contact
from .study import inputs,refinement_metrics,write


def check(path):
    path=ROOT/path;base=path.parent;s=json.loads(path.read_text())
    assert s['audit']==audit()
    verify_hashes(ROOT,s['source_sha256'])
    verify_hashes(base/'source_snapshot',s['source_sha256'])
    verify_hashes(base,s['artifacts_sha256'])
    old=frozen();D=old['batch_fit']['diffusivity_m2_s'];n=s['selection']['selected_sections']
    records={};selected=[]
    for key in s['case_ids']:
        r=json.loads((base/'cases'/f'{key}.json').read_text());p=r['inputs'];actual=contact(**p)
        for field,value in r['result'].items():
            if field=='total_s':continue
            if value is None:assert actual[field] is None
            else:close(actual[field],value,f'{key}: {field}')
        inventories({**r,'inputs':{**p,'active_inventory_fraction':1.}})
        assert r['metrics']==compare(actual,r['series'])
        assert p['diffusivity_m2_s']==D
        assert p['particle_density_kg_m3']==740/(1-.427)
        assert p['capacity_kg_kg']==.1115 and p['affinity_m3_kg']==139.
        assert 'active_inventory_fraction' not in p and 'axial_dispersion_m2_s' not in p
        assert r['result']['sections']==r['result']['axial_cells']==p['sections']
        if r['kind']!='contact_scenario':assert p['sections']==n
        expected_film=(2+1.1*(1000*p['superficial_velocity_m_s']*p['particle_diameter_m']/(1/1200))**.6*
            ((1/1200)/(1000*1e-9))**(1/3))*1e-9/p['particle_diameter_m']
        close(p['film_m_s'],expected_film,'Unfitted film correlation')
        Q=p['area_m2']*p['superficial_velocity_m_s']
        close(actual['liquid_residence_time_s'],p['porosity']*p['length_m']/p['superficial_velocity_m_s'],'Liquid residence time')
        close(actual['section_liquid_volume_m3']*n if p['sections']==n else actual['section_liquid_volume_m3']*p['sections'],
            p['porosity']*p['area_m2']*p['length_m'],'Section liquid volume sum')
        close(actual['section_dry_mass_kg']*p['sections'],actual['physical_adsorbent_mass_kg'],'Section mass sum')
        # Independently inspect SI CSV exports, including the outlet orientation.
        table=np.genfromtxt(base/'cases'/f'{key}.csv',delimiter=',',names=True)
        for name in table.dtype.names:close(table[name],actual[name],f'{key}: CSV {name}')
        records[key]=r
        if r['kind']=='selected':selected.append(r)
    # Every calibration objective is replayed using calibration observations only.
    t,y=observations('calibration');scores={}
    for i,row in enumerate(s['selection']['evaluations']):
        p=inputs(.0008,D,row['sections'],radial_cells=row['radial_cells'],dt_s=row['dt_s'])
        r=contact(**p)
        pred=np.interp(t,r['time_s'],np.array(r['outlet_concentration_kg_m3'])/.15)
        score=float(np.sqrt(np.mean((pred-y)**2)))
        close(score,row['RMSE'],'Replayed calibration score')
        close(r['normalized_mass_residual'],row['normalized_mass_residual'],'Search conservation')
        scores[(row['sections'],row['radial_cells'],row['dt_s'])]=score
    assert {key[0] for key in scores if key[1:]==(16,2.)}==set(range(1,65))
    fine={key[0]:value for key,value in scores.items() if key[1:]==(32,.5)}
    assert min(fine,key=fine.get)==n
    assert s['selection']['holdout_used'] is False
    assert s['selection']['fitted_curve']=='Figure 7, 0.8 mm only'
    refinement={};qualified_refinement={};qualified_t05={}
    for series in ('calibration','holdout'):
        a=records[f'selected_{series}']['result'];refinement[series]={};qualified_refinement[series]={}
        for label in ('radial','time'):
            refinement[series][label]=refinement_metrics(a,records[f'{label}_{series}']['result'])
            qualified_refinement[series]['late_'+label]=refinement_metrics(a,records[f'{label}_{series}']['result'],thresholds=('t50',),time_min_s=5.)
        early=records[f'early_{series}']['result'];qualified_t05[series]=records[f'early_{series}']['metrics']['t05']
        for label in ('radial','time'):
            qualified_refinement[series]['early_'+label]=refinement_metrics(early,records[f'early_{label}_{series}']['result'],thresholds=('t05',))
    assert refinement==s['refinement']
    assert qualified_refinement==s['qualified_refinement']
    assert qualified_t05==s['qualified_t05']
    numerical=all(v['passed'] for d in qualified_refinement.values() for v in d.values())
    assert numerical==s['numerical_refinement_passed']
    physical=all(r['metrics']['RMSE_screen_passed'] and r['metrics']['t50_screen_passed'] and
        r['metrics']['t05_screen_passed'] is not False for r in selected)
    assert physical==s['physical_status']['common_column_comparison_passed']
    assert s['physical_status']['batch_component_preserved']
    assert not s['physical_status']['process_energy_recommendation_supported']
    assert not s['physical_status']['experimental_pressure_validation']
    assert s['new_OpenLB_runs']==s['new_DWSIM_runs']==0
    assert len(selected)==2
    # Compact documentation fixture retains input snapshots and actual endpoints.
    compact={k:s[k] for k in ['run_id','physical_status','numerical_refinement_passed',
        'selection','refinement','qualified_refinement','qualified_t05','audit',
        'maximum_normalized_mass_residual','maximum_normalized_interface_flux_mismatch']}
    compact['cases']=[{k:r[k] for k in ['id','inputs','metrics','kind','series']} |
        {'result_scalars':{k:v for k,v in r['result'].items() if not isinstance(v,list) and k!='total_s'},
         'last_outlet_kg_m3':r['result']['outlet_concentration_kg_m3'][-1],
         'last_loading_kg_kg':r['result']['mean_loading_kg_kg'][-1]}
         for r in records.values() if r['kind'] in ('selected','packing_scenario')]
    write(ROOT/'docs/column_contact_reference.json',compact)
    for name in ('contact_comparison','inventory_audit'):
        shutil.copy2(base/f'{name}.png',ROOT/'docs'/f'column_{name}.png')
    message=dict(verified=True,forward_case_replays=len(records),
        calibration_evaluation_replays=len(scores),artifacts=len(s['artifacts_sha256']),
        source_files=len(s['source_sha256']),numerical_refinement_passed=numerical,
        physical_comparison_passed=physical)
    write(base/'verification.json',message)
    print(json.dumps(message,indent=2),flush=True)
    return message


def main():
    p=json.loads((ROOT/'results/adsorption/column_contact/latest.json').read_text())
    check(p['summary'])


if __name__=='__main__':main()
