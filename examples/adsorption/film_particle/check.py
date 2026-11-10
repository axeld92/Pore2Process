"""Independent inventories, source contracts and real Julia replay of each case."""
from pore2process_demo.provenance import verified_source_path
import json
import math
from pathlib import Path

import numpy as np

from examples.adsorption.alhamed.source import sha
from examples.adsorption.alhamed.study import compare
from .bridge import batch,bed,call
from .data import ROOT,evidence,observations
from .metrics import batch_metrics,physical_status
from .run import crossing_values,write


def verify_hashes(base,manifest):
    base=Path(base).resolve()
    if base == ROOT.resolve():
        for name,digest in manifest.items():
            verified_source_path(name,digest)
        return
    for name,digest in manifest.items():
        p=(base/name).resolve()
        if not p.is_relative_to(base) or sha(p)!=digest:
            raise ValueError(f'Changed or invalid artifact: {name}')


def close(actual,expected,label):
    np.testing.assert_allclose(actual,expected,rtol=1e-9,atol=1e-13,err_msg=label)


def inventories(record):
    p=record['inputs'];r=record['result'];is_batch=record['id'].startswith('batch_')
    n=r['radial_cells'];w=np.diff(np.linspace(0,1,n+1)**3)
    close(r['radial_volume_weights'],w,'Spherical shell weights')
    q=np.asarray(r['final_shell_loading_kg_kg']);c=np.asarray(r['final_concentration_kg_m3'])
    assert q.shape==(r['axial_cells'],n)
    assert np.all(q>=0) and np.all(c>=0)
    assert np.max(q)<=p['capacity_kg_kg']*(1+1e-10)
    mean=float(np.mean(q@w))
    if is_batch:
        physical=p['adsorbent_mass_kg'];active=physical
        fluid=p['fluid_volume_m3']-len(p['sampling_times_s'])*p['sample_volume_m3']
        close(r['fluid_volume_m3'][-1],fluid,'Bath sample-volume balance')
    else:
        physical=(1-p['porosity'])*p['particle_density_kg_m3']*p['area_m2']*p['length_m']
        active=physical*p['active_inventory_fraction']
        fluid=p['porosity']*p['area_m2']*p['length_m']
    close(r['physical_adsorbent_mass_kg'],physical,'Physical dry mass')
    close(r['active_adsorbent_mass_kg'],active,'Active dry mass')
    close(r['mean_loading_kg_kg'][-1],mean,'Particle average from shells')
    close(r['adsorbed_inventory_kg'][-1],active*mean,'Adsorbed inventory')
    close(r['liquid_inventory_kg'][-1],fluid*np.mean(c),'Liquid inventory')
    close(r['inventory_kg'][-1],active*mean+fluid*np.mean(c),'Final total inventory')
    close(r['adsorbed_inventory_kg'],active*np.array(r['mean_loading_kg_kg']),'Loading inventory history')
    residual=np.array(r['inventory_kg'])-r['inventory_kg'][0]-np.array(r['cumulative_inlet_kg'])+np.array(r['cumulative_outlet_kg'])+np.array(r['cumulative_sampled_solute_kg'])
    close(r['mass_residual_kg'],residual,'Reconstructed mass residual')
    normalized=np.max(np.abs(residual))/max(r['inventory_kg'][0],r['cumulative_inlet_kg'][-1],1e-30)
    close(normalized,r['normalized_mass_residual'],'Normalized conservation')
    assert normalized<=1e-8 and r['normalized_interface_flux_mismatch']<=1e-8
    if is_batch:
        assert r['cumulative_inlet_kg'][-1]==0 and r['cumulative_outlet_kg'][-1]==0
        removed=0.
        times=np.array(r['time_s'])
        for t in p['sampling_times_s']:
            i=int(np.argmin(abs(times-t)));assert abs(times[i]-t)<=1e-8*max(t,1)
            removed+=p['sample_volume_m3']*r['outlet_concentration_kg_m3'][i]
        close(removed,r['cumulative_sampled_solute_kg'][-1],'Serial sample solute loss')
    else:
        dt=np.diff(r['time_s']);Q=p['superficial_velocity_m_s']*p['area_m2']
        expected_in=np.r_[0,np.cumsum(Q*p['inlet_concentration_kg_m3']*dt)]
        expected_out=np.r_[0,np.cumsum(Q*np.asarray(r['outlet_concentration_kg_m3'])[1:]*dt)]
        close(r['cumulative_inlet_kg'],expected_in,'Integrated feed')
        close(r['cumulative_outlet_kg'],expected_out,'Integrated effluent')
    assert p['capacity_kg_kg']==.1115 and p['affinity_m3_kg']==139.


def check(path):
    path=Path(path);base=path.parent;s=json.loads(path.read_text())
    assert s['evidence']==evidence()
    verify_hashes(ROOT,s['source_sha256']);verify_hashes(base/'source_snapshot',s['source_sha256'])
    verify_hashes(base,s['artifacts_sha256'])
    obs=observations();records={}
    for key in s['case_ids']:
        record=json.loads((base/'cases'/f'{key}.json').read_text());records[key]=record
        actual=(batch if key.startswith('batch_') else bed)(**record['inputs'])
        for field,value in record['result'].items():
            if field=='total_s':continue
            close(actual[field],value,f'{key}: {field}')
        inventories(record)
        if 'metrics' in record:
            if key.startswith('batch_'):
                label=record['metadata'].get('diameter_id','d1500')
                m=batch_metrics(actual,obs[label])
            else:
                m=compare(actual,record['metadata']['series'])
            assert m==record['metrics']
    endpoint=json.loads((base/'endpoint_audit.json').read_text())
    for row in endpoint['source_endpoints']:
        actual=call('bath_equilibrium',**row['inputs'])
        for key,value in row['result'].items():close(actual[key],value,'Finite-bath endpoint replay')
        p=row['inputs'];c0=p['initial_concentration_kg_m3'];beta=p['adsorbent_mass_kg']/p['fluid_volume_m3']
        b=p['affinity_m3_kg'];qm=p['capacity_kg_kg'];B=1+beta*qm*b-c0*b
        analytical=2*c0/(B+math.sqrt(B*B+4*b*c0))
        close(actual['concentration_kg_m3'],analytical,'Independent quadratic endpoint')
    for row in s['initial_rate_audit']:
        value=call('initial_film_rate',**row['inputs'])['initial_rate_kg_kg_s']
        close(value,row['predicted_rate_kg_kg_s'],'Initial film-supply limit')
        close(value*60000,row['predicted_rate_mg_g_min'],'Initial-rate units')
    batches=[records[f'batch_{key}'] for key in obs]
    columns=[records[key] for key in s['case_ids'] if records[key]['kind']=='column']
    assert physical_status(endpoint,batches,columns)==s['physical_status']
    fit=s['batch_fit'];anchor=s['column_film_fit']
    assert set(fit['fitted_diameter_ids'])=={'d1500','d800'}
    assert anchor['fitted_curve']=='Figure 7, 0.8 mm only'
    for record in batches:
        p=record['inputs']
        close(p['diffusivity_m2_s'],fit['diffusivity_m2_s'],'Shared batch Ds')
        close(p['film_m_s'],fit['batch_film_ref_m_s']*.0008/p['particle_diameter_m'],'Shared film size law')
    for record in columns:
        p=record['inputs'];meta=record['metadata']
        close(p['diffusivity_m2_s'],fit['diffusivity_m2_s'],'Frozen column Ds')
        assert p['particle_density_kg_m3']==740/(1-.427)
        velocity=p['superficial_velocity_m_s']/(p['porosity'] if meta['velocity_convention']=='interstitial' else 1)
        Re=1000*velocity*p['particle_diameter_m']/(1/1200);Sc=(1/1200)/(1000*1e-9)
        expected=meta['multiplier']*(2+1.1*Re**.6*Sc**(1/3))*1e-9/p['particle_diameter_m']
        close(p['film_m_s'],expected,'Independent column-film correlation')
        if meta['family']=='anchored':close(meta['multiplier'],anchor['multiplier'],'One common film multiplier')
    qnorm=endpoint['source_endpoints'][0]['result']['loading_kg_kg']
    numerical=True
    for key,data in obs.items():
        r=records[f'batch_{key}']['result'];baseline=np.interp(data['time_s'],r['time_s'],r['mean_loading_kg_kg'])
        for name in ['radial','time']:
            r=records[f'batch_{key}_{name}']['result'];q=np.interp(data['time_s'],r['time_s'],r['mean_loading_kg_kg'])
            change=float(np.max(abs(q-baseline))/qnorm)
            close(change,s['batch_refinement'][key][name],'Batch refinement')
            numerical &= change<=.02
    for key,changes in s['column_refinement'].items():
        baseline=crossing_values(records[key]['result'])
        for name,expected in changes.items():
            values=crossing_values(records[f'{key}_{name}']['result'])
            for metric in ['t05','t50']:
                value=None if baseline[metric] is None or values[metric] is None else abs(values[metric]/baseline[metric]-1)
                if value is None:assert expected[metric] is None
                else:close(value,expected[metric],'Column refinement')
                if metric=='t50' or key=='primary_holdout':
                    numerical &= value is not None and value<=.02
    for family,expected in s['early_calibration_refinement'].items():
        values={key:crossing_values(records[f'{family}_early_{key}']['result'])['t05']
                for key in ['base','axial','radial','time']}
        for key,value in values.items():close(value,expected['crossings_s'][key],'Short-window crossing')
        for key in ['axial','radial','time']:
            change=abs(values[key]/values['base']-1)
            close(change,expected['relative_changes'][key],'Short-window refinement')
            numerical &= change<=.02
    assert bool(numerical)==s['numerical_refinement_passed']
    return dict(status='passed_source_inventory_and_replay_checks',Julia_cases_replayed=len(records),
        independent_endpoint_checks=len(endpoint['source_endpoints']),artifacts_checked=len(s['artifacts_sha256']),
        source_files_checked=len(s['source_sha256']),numerical_refinement_passed=bool(numerical),
        physical_status=s['physical_status'],new_OpenLB_runs=0,new_DWSIM_runs=0)


def main():
    pointer=ROOT/'results/adsorption/film_particle/latest.json'
    path=Path(json.loads(pointer.read_text())['summary']);result=check(path)
    write(path.parent/'verification.json',result)
    s=json.loads(path.read_text());reference={k:s[k] for k in ['run_id','versions','batch_fit','column_film_fit',
        'physical_status','numerical_refinement_passed','batch_refinement','column_refinement',
        'early_calibration_refinement','endpoint_audit','initial_rate_audit','batch_metrics','column_metrics','limitations']}
    reference['summary_path']=str(path);reference['verification']=result
    chosen=[f'batch_{key}' for key in observations()]+['primary_calibration','primary_holdout','anchored_calibration','anchored_holdout']
    reference['cases']=[]
    for key in chosen:
        case=json.loads((path.parent/'cases'/f'{key}.json').read_text())
        r=case['result'];case['result_scalars']={k:v for k,v in r.items() if not isinstance(v,list) and k!='total_s'}
        case['last_loading_kg_kg']=r['mean_loading_kg_kg'][-1]
        case['last_outlet_kg_m3']=r['outlet_concentration_kg_m3'][-1]
        del case['result'];reference['cases'].append(case)
    write(ROOT/'examples/adsorption/film_particle/reference/film_particle_reference.json',reference)
    for name in ['batch_comparison.png','column_comparison.png']:
        import shutil
        shutil.copy2(path.parent/name,ROOT/'examples/adsorption/film_particle/figures'/f'film_particle_{name}')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
