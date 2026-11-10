"""One contact topology chosen on 0.8 mm only, then frozen for 0.45 mm."""
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
from time import perf_counter
from uuid import uuid4
import numpy as np

from adsorption.alhamed.source import ROOT,sha
from adsorption.alhamed.study import observations,compare,crossing
from adsorption.film_particle.cases import column_inputs
from pore2process_demo.bridge import runtime_versions
from .audit import audit,FROZEN_RUN,PROTOCOL
from .bridge import module,contact


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def inputs(diameter,D,sections,*,porosity=.427,radial_cells=32,dt_s=.5):
    p,film=column_inputs(diameter,D,radial_cells=radial_cells,dt_s=dt_s)
    for key in ('active_inventory_fraction','axial_dispersion_m2_s','axial_cells'):
        del p[key]
    p.update(sections=int(sections),porosity=float(porosity))
    return p


def select_sections(D):
    # Holdout is never read here. Exactly one discrete choice is selected.
    t,y=observations('calibration');evaluations=[];cache={}
    def evaluate(n,nr,dt):
        key=(n,nr,dt)
        if key not in cache:
            p=inputs(.0008,D,n,radial_cells=nr,dt_s=dt);r=contact(**p)
            pred=np.interp(t,r['time_s'],np.asarray(r['outlet_concentration_kg_m3'])/.15)
            score=float(np.sqrt(np.mean((pred-y)**2)))
            if r['normalized_mass_residual']>1e-8:raise ValueError('Search conservation failed')
            cache[key]=score
            evaluations.append(dict(sections=n,radial_cells=nr,dt_s=dt,RMSE=score,
                normalized_mass_residual=r['normalized_mass_residual']))
        return cache[key]
    selected=min(range(1,65),key=lambda n:evaluate(n,16,2.))
    print('Coarse selected sections',selected,flush=True)
    nearby=range(max(1,selected-1),min(64,selected+1)+1)
    fine=min(nearby,key=lambda n:evaluate(n,32,.5))
    if fine!=selected:
        fine=min(range(1,65),key=lambda n:evaluate(n,32,.5))
    return dict(selected_sections=int(fine),fitted_curve='Figure 7, 0.8 mm only',
        holdout_used=False,range=[1,64],boundary_selection=fine in (1,64),
        fine_candidate_RMSE=cache[(fine,32,.5)],evaluations=evaluations,
        changed_coarse_selection=fine!=selected,
        particle_film_equilibrium_inventory_fitted=False)


def save_case(base,name,p,series,kind):
    r=contact(**p)
    if r['normalized_mass_residual']>1e-8 or r['normalized_interface_flux_mismatch']>1e-8:
        raise ValueError(f'{name}: conservation/interface failed')
    record=dict(id=name,inputs=p,result=r,series=series,kind=kind,metrics=compare(r,series))
    write(base/'cases'/f'{name}.json',record)
    with (base/'cases'/f'{name}.csv').open('w',newline='') as stream:
        names=['time_s','outlet_concentration_kg_m3','mean_loading_kg_kg',
            'inventory_kg','adsorbed_inventory_kg','liquid_inventory_kg',
            'cumulative_inlet_kg','cumulative_outlet_kg','mass_residual_kg']
        writer=csv.writer(stream);writer.writerow(names)
        writer.writerows(zip(*(r[k] for k in names)))
    return record


def refinement_metrics(a,b,*,thresholds=('t05','t50'),time_min_s=0.):
    out={}
    for label,threshold in [('t05',.05),('t50',.5)]:
        if label not in thresholds:continue
        x=crossing(a['time_s'],np.array(a['outlet_concentration_kg_m3'])/.15,threshold)
        y=crossing(b['time_s'],np.array(b['outlet_concentration_kg_m3'])/.15,threshold)
        out[label+'_relative_change']=None if x is None or y is None else abs(y/x-1)
    t=np.array(a['time_s']);mask=t>=time_min_s
    out['maximum_fraction_change']=float(np.max(abs(np.interp(t[mask],b['time_s'],b['outlet_concentration_kg_m3'])-
                                                      np.array(a['outlet_concentration_kg_m3'])[mask])))/.15
    out['time_window_s']=[float(time_min_s),float(a['time_s'][-1])]
    out['passed']=all(out[k+'_relative_change'] is not None and out[k+'_relative_change']<=.02 for k in thresholds) and out['maximum_fraction_change']<=.005
    return out


def figures(base,records,audit_data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained',sharey=True)
    for ax,series,diam in zip(axes,('calibration','holdout'),(.8,.45)):
        t,y=observations(series);ax.plot(t/60,y,'o',ms=4,mfc='none',color='black',label='Paper markers')
        for r in records:
            if r['series']!=series or r['kind'] not in ('selected','contact_scenario'):continue
            n=r['inputs']['sections'];selected=r['kind']=='selected'
            ax.plot(np.array(r['result']['time_s'])/60,np.array(r['result']['outlet_concentration_kg_m3'])/.15,
                lw=2.2 if selected else 1,ls='-' if selected else '--',alpha=1 if selected else .6,
                label=f'Selected: {n} section(s)' if selected else f'Scenario: {n} sections')
        ax.set(xlim=(0,300),ylim=(0,1.02),xlabel='Time (min)',
            title=f'{diam:g} mm: '+('calibration' if series=='calibration' else 'withheld'))
        ax.grid(alpha=.2);ax.legend(loc='lower right',fontsize=8)
    axes[0].set_ylabel('Outlet / inlet phenol concentration (1)')
    fig.suptitle('Common mixed-section model; batch particles and printed inventory fixed')
    fig.savefig(base/'contact_comparison.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    ax=axes[0]
    for series in ('calibration','holdout'):
        t,y=observations(series);tt=np.r_[0,t];deficit=1-np.r_[0,y]
        retained=np.r_[0,np.cumsum(np.diff(tt)*(deficit[:-1]+deficit[1:])/2)]*7.3e-7*.15
        ax.plot(tt/60,retained*1e6,label='.8 mm' if series=='calibration' else '.45 mm')
    ax.axhline(audit_data['equilibrium']['full_equilibrium_inventory_kg']*1e6,color='black',ls='--',label='Full equilibrium target')
    ax.set(xlabel='Time (min)',ylabel='Retained phenol (mg)',title='Digitized window leaves the tail unresolved')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1]
    for series in ('calibration','holdout'):
        values=[r for r in records if r['series']==series and r['kind'] in ('selected','packing_scenario')]
        values.sort(key=lambda r:r['inputs']['porosity'])
        ax.plot([r['inputs']['porosity'] for r in values],
            [r['metrics']['t50']['predicted_s']/60 for r in values],'o-',label='.8 mm' if series=='calibration' else '.45 mm')
    ax.set(xlabel='Chosen bed voidage (1)',ylabel='Time to 50% outlet concentration (min)',title='Packing scenarios, no refit')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(base/'inventory_audit.png',dpi=180);plt.close(fig)


def main():
    start=perf_counter();a=audit();cold=perf_counter();module();cold=perf_counter()-cold
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=ROOT/'results/adsorption/column_contact/runs'/run_id;(base/'cases').mkdir(parents=True)
    sources=[*Path(__file__).parent.glob('*.py'),PROTOCOL,FROZEN_RUN,
        ROOT/'julia/adsorption/ColumnContact.jl',ROOT/'julia/adsorption/test_column_contact.jl',
        ROOT/'tests/test_column_contact.py',
        *[ROOT/name for name in a['prior_batch_source_sha256']]]
    sources={str(p.relative_to(ROOT)):sha(p) for p in sources}
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    write(base/'audit.json',a);write(base/'source_sha256.json',sources)
    D=a['batch_diffusivity_m2_s'];fit=select_sections(D);write(base/'contact_selection.json',fit)
    n=fit['selected_sections'];records=[];refinement={};qualified_refinement={};qualified_t05={}
    for series,d in [('calibration',.0008),('holdout',.00045)]:
        p=inputs(d,D,n);r=save_case(base,f'selected_{series}',p,series,'selected');records.append(r)
        print(series,r['metrics'],flush=True)
        refinement[series]={}
        qualified_refinement[series]={}
        for label,update in [('radial',{'radial_cells':64}),('time',{'dt_s':.25})]:
            rr=save_case(base,f'{label}_{series}',{**p,**update},series,'refinement');records.append(rr)
            refinement[series][label]=refinement_metrics(r['result'],rr['result'])
            qualified_refinement[series]['late_'+label]=refinement_metrics(r['result'],rr['result'],thresholds=('t50',),time_min_s=5.)
        # Cover the omitted startup with a separately resolved 0..5 s window.
        # Preserve the failed coarse t05 checks; do not relax their tolerances.
        early_p={**p,'t_end_s':5.,'radial_cells':64,'dt_s':.002}
        early=save_case(base,f'early_{series}',early_p,series,'early_refinement');records.append(early)
        qualified_t05[series]=early['metrics']['t05']
        for label,update in [('radial',{'radial_cells':128}),('time',{'dt_s':.001})]:
            rr=save_case(base,f'early_{label}_{series}',{**early_p,**update},series,'early_refinement');records.append(rr)
            qualified_refinement[series]['early_'+label]=refinement_metrics(early['result'],rr['result'],thresholds=('t05',))
        for scenario_n in (1,4,16,64):
            if scenario_n==n:continue
            records.append(save_case(base,f'sections{scenario_n}_{series}',inputs(d,D,scenario_n),series,'contact_scenario'))
        for eps in (.35,.5,.6):
            records.append(save_case(base,f'porosity{eps}_{series}',inputs(d,D,n,porosity=eps),series,'packing_scenario'))
    figures(base,records,a)
    physical=all(r['metrics']['RMSE_screen_passed'] and r['metrics']['t50_screen_passed'] and
        r['metrics']['t05_screen_passed'] is not False for r in records if r['kind']=='selected')
    summary=dict(schema='pore2process.column-contact.v1',run_id=run_id,timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sources,audit=a,selection=fit,case_ids=[r['id'] for r in records],
        refinement=refinement,qualified_refinement=qualified_refinement,qualified_t05=qualified_t05,
        numerical_refinement_passed=all(v['passed'] for s in qualified_refinement.values() for v in s.values()),
        physical_status=dict(batch_component_preserved=True,common_column_comparison_passed=physical,
            experimental_pressure_validation=False,process_energy_recommendation_supported=False),
        model='N equal mixed sections; full active solid and liquid storage; zero inter-section dispersion',
        versions=runtime_versions(),cold_initialization_s=cold,study_total_s=perf_counter()-start,
        maximum_normalized_mass_residual=max(r['result']['normalized_mass_residual'] for r in records),
        maximum_normalized_interface_flux_mismatch=max(r['result']['normalized_interface_flux_mismatch'] for r in records),
        new_OpenLB_runs=0,new_DWSIM_runs=0)
    summary['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write(base/'summary.json',summary)
    write(base.parent.parent/'latest.json',{'summary':str((base/'summary.json').relative_to(ROOT))})
    print('Saved',base,'numerical refinement',summary['numerical_refinement_passed'],'physical comparison',physical,flush=True)


if __name__=='__main__':main()
