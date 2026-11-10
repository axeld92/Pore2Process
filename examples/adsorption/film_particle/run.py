"""Reproducible finite-bath calibration and withheld column prediction."""
from pore2process_demo.provenance import implementation_files, source_manifest
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.bridge import runtime_versions
from examples.adsorption.alhamed.source import sha
from examples.adsorption.alhamed.study import compare, crossing
from .bridge import module,batch,bed,call
from .cases import batch_inputs,column_inputs
from .data import ROOT,PIXELS,PROTOCOL,evidence,observations
from .fitting import endpoint_audit,fit_batch,fit_column_film
from .metrics import batch_metrics,physical_status
from .plot import figures


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def save_case(base,name,function,inputs,kind,metrics=None,metadata=None):
    result=function(**inputs)
    if result['normalized_mass_residual']>1e-8 or result['normalized_interface_flux_mismatch']>1e-8:
        raise ValueError(f'{name}: numerical conservation or interface screen failed.')
    record=dict(id=name,kind=kind,inputs=inputs,result=result,metadata=metadata or {})
    if metrics is not None:
        record['metrics']=metrics(result)
    write(base/'cases'/f'{name}.json',record)
    with (base/'cases'/f'{name}_series.csv').open('w',newline='') as stream:
        names=['time_s','outlet_concentration_kg_m3','mean_loading_kg_kg','inventory_kg',
            'adsorbed_inventory_kg','liquid_inventory_kg','cumulative_inlet_kg',
            'cumulative_outlet_kg','cumulative_sampled_solute_kg','mass_residual_kg']
        writer=csv.writer(stream);writer.writerow(names)
        writer.writerows(zip(*(result[k] for k in names)))
    return record


def crossing_values(r):
    y=np.array(r['outlet_concentration_kg_m3'])/.15
    return {k:crossing(r['time_s'],y,p) for k,p in [('t05',.05),('t50',.5)]}


def main():
    start=perf_counter();module();initialization=perf_counter()-start
    e=evidence();obs=observations()
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=ROOT/'results/adsorption/film_particle/runs'/run_id
    (base/'cases').mkdir(parents=True)
    source_files=[*Path(__file__).parent.glob('*.py'),ROOT/'julia/adsorption/FilmParticle.jl',
        ROOT/'julia/adsorption/test_film_particle.jl',PIXELS,PROTOCOL,
        ROOT/'examples/adsorption/alhamed/source.py',ROOT/'examples/adsorption/alhamed/study.py',
        ROOT/'literature/sources/alhamed2009/figure7_pixels.csv',ROOT/'src/pore2process_demo/bridge.py',
        ROOT/'uv.lock',ROOT/'julia/Project.toml',ROOT/'julia/Manifest.toml']
    source_files = implementation_files(source_files)
    sources=source_manifest(source_files)
    for p in source_files:
        target=base/'source_snapshot'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    write(base/'evidence.json',e)
    print('Auditing finite-bath equilibrium...',flush=True)
    endpoint=endpoint_audit();write(base/'endpoint_audit.json',endpoint)
    print('Fitting two shared batch parameters, with small particles withheld...',flush=True)
    fit=fit_batch();write(base/'batch_fit.json',fit)
    D=fit['diffusivity_m2_s'];kf=fit['batch_film_ref_m_s']
    rates=[]
    for row in e['batch_table_3']:
        diameter=row['diameter_mm']*.001
        p=dict(particle_diameter_m=diameter,particle_density_kg_m3=740/(1-.427),
            film_m_s=kf*.0008/diameter,concentration_kg_m3=row['inlet_mg_L']*.001)
        value=call('initial_film_rate',**p)['initial_rate_kg_kg_s']
        rates.append(dict(experiment=row['experiment'],inputs=p,predicted_rate_kg_kg_s=value,
            predicted_rate_mg_g_min=value*60000,reported_PSO_fit_h_mg_g_min=row['fitted_initial_rate_mg_g_min'],
            role='descriptive comparison to a published extrapolated fit; not a measured t=0 rate or a fitting target'))
    write(base/'initial_rate_audit.json',rates)
    print('Batch fit',D,kf,'search RMSE mg/g',fit['search_RMSE_mg_g'],flush=True)
    records=[];batch_records=[];column_records=[];refinement={};batch_refinement={}
    for key,data in obs.items():
        p=batch_inputs(data['diameter_m'],D,kf)
        record=save_case(base,f'batch_{key}',batch,p,'batch',lambda r,d=data:batch_metrics(r,d),{'diameter_id':key,'role':data['role']})
        records.append(record);batch_records.append(record)
        q=np.interp(data['time_s'],record['result']['time_s'],record['result']['mean_loading_kg_kg'])
        changes={}
        for name,update in [('radial',{'radial_cells':48}),('time',{'dt_s':15.})]:
            refined=save_case(base,f'batch_{key}_{name}',batch,{**p,**update},'batch_refinement')
            records.append(refined)
            qr=np.interp(data['time_s'],refined['result']['time_s'],refined['result']['mean_loading_kg_kg'])
            source_loading=endpoint['source_endpoints'][0]['result']['loading_kg_kg']
            changes[name]=float(np.max(np.abs(qr-q))/source_loading)
        batch_refinement[key]=changes
        sample=save_case(base,f'batch_{key}_sampling',batch,
            {**p,'sampling_times_s':data['time_s'].tolist(),'sample_volume_m3':1e-6},'batch_sampling')
        records.append(sample)
        print(key,'RMSE mg/g',record['metrics']['RMSE_mg_g'],flush=True)
    label=save_case(base,'batch_d1470_label_scenario',batch,batch_inputs(.00147,D,kf),'diameter_label_scenario',
        lambda r:batch_metrics(r,obs['d1500']),{'printed_labels_mm':[1.47,1.5]})
    records.append(label)
    print('One-anchor column film fit; no inventory, equilibrium or diffusivity refit...',flush=True)
    anchor=fit_column_film(D);write(base/'column_film_fit.json',anchor)
    families=[('primary',{}),('interstitial',{'interstitial':True}),('dispersion',{'dispersion_factor':.5}),
        ('anchored',{'multiplier':anchor['multiplier']}),
        ('active50',{'multiplier':anchor['multiplier'],'active_fraction':.5}),
        ('active75',{'multiplier':anchor['multiplier'],'active_fraction':.75})]
    for family,options in families:
        for series,diameter in [('calibration',.0008),('holdout',.00045)]:
            p,film=column_inputs(diameter,D,**options)
            r=save_case(base,f'{family}_{series}',bed,p,'column',lambda r,s=series:compare(r,s),
                {**film,'family':family,'series':series,'column_fit_used':family not in ['primary','interstitial','dispersion']})
            records.append(r);column_records.append(r)
            if family in ['primary','anchored']:
                changes={};baseline=crossing_values(r['result'])
                for name,update in [('radial',{'radial_cells':32}),('axial',{'axial_cells':160}),('time',{'dt_s':.5})]:
                    refined=save_case(base,f'{family}_{series}_{name}',bed,{**p,**update},'column_refinement')
                    records.append(refined);values=crossing_values(refined['result'])
                    changes[name]={k:None if baseline[k] is None or values[k] is None else abs(values[k]/baseline[k]-1)
                                   for k in ['t05','t50']}
                refinement[f'{family}_{series}']=changes
            print(family,series,'RMSE',r['metrics']['RMSE'],'t50 min',r['metrics']['t50']['predicted_s']/60,flush=True)
    # Independently resolve the sub-minute calibration t05 on a short window;
    # retaining full-run coarse values makes the distinction explicit.
    early={}
    early_specs=[('primary',1.,.0008,320,640,24,.025,120.),
        ('anchored',anchor['multiplier'],.0008,320,640,24,.025,120.),
        ('anchored_holdout',anchor['multiplier'],.00045,1000,2000,32,.25,1500.)]
    for family,multiplier,diameter,nz,nz_ref,nr,dt,tend in early_specs:
        local={};early_records={}
        for name,n_cells,r_cells,h in [('base',nz,nr,dt),('axial',nz_ref,nr,dt),
                                      ('radial',nz,2*nr,dt),('time',nz,nr,dt/2)]:
            p,_=column_inputs(diameter,D,multiplier=multiplier,axial_cells=n_cells,radial_cells=r_cells,dt_s=h)
            p['t_end_s']=tend
            r=save_case(base,f'{family}_early_{name}',bed,p,'early_column_refinement')
            records.append(r);value=crossing_values(r['result'])['t05']
            early_records[name]=value
        local['crossings_s']=early_records
        local['relative_changes']={k:abs(v/early_records['base']-1) for k,v in early_records.items() if k!='base'}
        early[family]=local
    # Near-optimal batch-film profiles are transfer scenarios, not confidence bounds.
    compatible=[x for x in fit['kf_profiles'] if x['RMSE_mg_g']<=fit['search_RMSE_mg_g']+.5]
    profile_cases=[]
    if compatible:
        selected={min(compatible,key=lambda x:x['best_Ds_m2_s'])['log10_batch_kf_ref'],
                  max(compatible,key=lambda x:x['best_Ds_m2_s'])['log10_batch_kf_ref']}
        for row in compatible:
            if row['log10_batch_kf_ref'] not in selected:continue
            for series,diameter in [('calibration',.0008),('holdout',.00045)]:
                p,film=column_inputs(diameter,row['best_Ds_m2_s'])
                r=save_case(base,f"profile_kf{row['log10_batch_kf_ref']}_{series}",bed,p,'batch_profile_transfer',
                    lambda r,s=series:compare(r,s),{'profile':row,'series':series,**film})
                records.append(r);profile_cases.append(r['id'])
    figures(base,batch_records,column_records,obs)
    numerical=all(v<=.02 for c in batch_refinement.values() for v in c.values())
    numerical &= all(v is not None and v<=.02 for key,c in refinement.items() for r in c.values()
                     for k,v in r.items() if k=='t50' or k=='t05' and key=='primary_holdout')
    numerical &= all(v<=.02 for c in early.values() for v in c['relative_changes'].values())
    science=physical_status(endpoint,batch_records,column_records)
    summary=dict(status='completed_comparison',run_id=run_id,evidence=e,source_sha256=sources,
        versions=runtime_versions(),Julia_initialization_s=initialization,elapsed_s=perf_counter()-start,
        batch_fit={k:v for k,v in fit.items() if k!='history'},
        column_film_fit={k:v for k,v in anchor.items() if k!='history'},
        endpoint_audit=endpoint,initial_rate_audit=rates,physical_status=science,numerical_refinement_passed=bool(numerical),
        batch_refinement=batch_refinement,column_refinement=refinement,early_calibration_refinement=early,
        refinement_scope='Full-run t50 and primary held-out t05; separate short windows for both calibration t05 and the sensitive anchored held-out t05.',
        batch_metrics={r['id']:r['metrics'] for r in batch_records},column_metrics={r['id']:r['metrics'] for r in column_records},
        profile_scenario_cases=profile_cases,case_ids=[r['id'] for r in records],Julia_final_cases=len(records),
        new_OpenLB_runs=0,new_DWSIM_runs=0,
        limitations=['Same-paper holdouts, not an independent experimental campaign.',
            'Homogeneous spherical loading diffusion neglects intraparticle liquid storage and real granule geometry.',
            'Fitted transport coefficients and the batch film-size scaling are effective assumptions.',
            'Source bulk density times dimensions is a derived inventory, not a measured per-run dry mass.',
            'Active fractions and dispersion choices are scenarios, not measured uncertainties.',
            'A passed numerical replay does not promote a failed column prediction.'])
    summary['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write(base/'summary.json',summary)
    write(ROOT/'results/adsorption/film_particle/latest.json',{'summary':str(base/'summary.json')})
    print(json.dumps({'summary':str(base/'summary.json'),'numerical_refinement_passed':bool(numerical),
                      'physical_status':science,'final_cases':len(records)},indent=2),flush=True)


if __name__=='__main__':
    main()
