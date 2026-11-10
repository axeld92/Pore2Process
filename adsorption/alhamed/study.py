"""Fit one published curve, then test diameter transfer with the real Julia solver."""
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import time
from uuid import uuid4

import numpy as np

from adsorption.bridge import column
from adsorption.bridge import module
from pore2process_demo.bridge import runtime_versions
from .source import ROOT, PDF, PIXELS, evidence, points, sha


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def crossing(t, y, threshold):
    t=np.asarray(t);y=np.asarray(y)
    found=np.flatnonzero(y>=threshold)
    if not len(found):
        return None
    i=int(found[0])
    if i==0:
        return float(t[0]) if y[0]==threshold else None
    return float(t[i-1]+(threshold-y[i-1])/(y[i]-y[i-1])*(t[i]-t[i-1]))


def observations(series, shift=(0.,0.)):
    p=points()[series]
    # Add a declared clean-bed origin only for interpolating crossings. It is
    # not an extra experimental point and is never included in the fit.
    return (np.maximum(0.,p['time_s']+shift[0]),
            np.clip(p['fraction']+shift[1],0.,1.))


def observed_crossing(series, threshold, shift=(0.,0.)):
    t,y=observations(series,shift)
    return crossing(np.r_[0.,t],np.r_[0.,y],threshold)


def inputs(rate, dispersion, cells=80, dt=.25):
    e=evidence();f=e['figure_7']
    return dict(length_m=f['length_m'],area_m2=f['area_m2'],porosity=.427,
        particle_density_kg_m3=740/(1-.427),
        superficial_velocity_m_s=f['flow_m3_s']/f['area_m2'],
        inlet_concentration_kg_m3=.150,capacity_kg_kg=.1115,affinity_m3_kg=139.,
        ldf_s_inv=float(rate),axial_dispersion_m2_s=float(dispersion),
        t_end_s=18000.,cells=cells,dt_s=dt,initial_concentration_kg_m3=0.,
        initial_loading_kg_kg=0.,linear_partition_m3_kg=None,equilibrium=False,
        breakthrough_fraction=.05)


def fit(shift=(0.,0.), start=None):
    """Bounded log-space coarse search and deterministic pattern refinement.

    Python only coordinates evaluations. Every forward prediction is produced
    by the existing Julia conservative finite-volume column model.
    """
    t,y=observations('calibration',shift)
    cache={};history=[]
    lower=np.array([-7.,-10.]);upper=np.array([-2.,-2.])

    def objective(x):
        x=np.clip(np.asarray(x),lower,upper)
        key=tuple(np.round(x,10))
        if key not in cache:
            # Search at 1 s, then independently evaluate/refine the final
            # predictions at 0.25 / 0.125 s, including the early crossing.
            params=inputs(10.**x[0],10.**x[1],dt=1.)
            result=column(**params)
            prediction=np.interp(t,result['time_s'],np.array(result['outlet_concentration_kg_m3'])/.15)
            error=float(np.mean((prediction-y)**2))
            cache[key]=error
            history.append(dict(log10_rate=float(x[0]),log10_dispersion=float(x[1]),
                MSE=error,normalized_mass_residual=result['normalized_mass_residual']))
        return cache[key]

    if start is None:
        seeds=[np.array([r,d]) for r in [-7,-6,-5,-4,-3,-2] for d in [-10,-8,-6,-4,-2]]
        x=min(seeds,key=objective)
    else:
        x=np.array(start,dtype=float)
    best=objective(x)
    directions=[np.array([a,b]) for a in [-1,0,1] for b in [-1,0,1] if a or b]
    for step in [1.,.5,.25,.125,.0625,.03125,.015625,.0078125]:
        for iteration in range(200):
            candidates=[np.clip(x+step*d,lower,upper) for d in directions]
            candidate=min(candidates,key=objective);score=objective(candidate)
            if score>=best-1e-14:
                break
            x=candidate;best=score
        else:
            raise RuntimeError(f'Pattern search did not settle: x={x}, MSE={best}; no fit accepted.')
        print(f'Search step {step:g}: log10 parameters {x}; RMSE {math.sqrt(best):.6g}',flush=True)
    return dict(log10_parameters=x.tolist(),effective_ldf_s_inv=float(10.**x[0]),
        dispersion_m2_s=float(10.**x[1]),RMSE=math.sqrt(best),
        fitted_series='calibration only; 0.8 mm',evaluations=len(history),
        search_cells=80,search_dt_s=1.,
        log10_bounds=[lower.tolist(),upper.tolist()],history=history,
        search_resolution_log10=.0078125,
        parameters_at_bounds=bool(np.any(np.minimum(x-lower,upper-x)<.008)))


def compare(result, series):
    t,y=observations(series)
    fraction=np.array(result['outlet_concentration_kg_m3'])/.15
    predicted=np.interp(t,result['time_s'],fraction)
    out=dict(RMSE=float(np.sqrt(np.mean((predicted-y)**2))),
             maximum_absolute_fraction_error=float(np.max(np.abs(predicted-y))),
             numerical_normalized_mass_residual=result['normalized_mass_residual'])
    e=evidence()['figure_7'];tx=e['reading_sensitivity_time_s'];dy=e['reading_sensitivity_fraction']
    for label,threshold in [('t05',.05),('t50',.5)]:
        obs=observed_crossing(series,threshold)
        pred=crossing(result['time_s'],fraction,threshold)
        corners=[observed_crossing(series,threshold,(sx*tx,sy*dy))
                 for sx in [-1,1] for sy in [-1,1]]
        valid=[x for x in corners if x is not None]
        bound=[min(valid),max(valid)] if valid else [None,None]
        out[label]=dict(observed_interpolated_s=obs,predicted_s=pred,
            reading_sensitivity_s=bound,
            relative_error=None if not obs or pred is None else abs(pred/obs-1),
            observed_crossing_resolved=bool(valid and min(valid)>0))
    out['RMSE_screen_passed']=out['RMSE']<=.1
    out['t50_screen_passed']=out['t50']['relative_error'] is not None and out['t50']['relative_error']<=.2
    c=out['t05']
    out['t05_screen_passed']=None if not c['observed_crossing_resolved'] else c['relative_error'] is not None and c['relative_error']<=.25
    return out


def source_equation_audit():
    """Independent evaluation of the printed fitted expression, not a solver.

    A square-root alternative is labelled explicitly, never silently inserted.
    No value from this audit calibrates the Julia forward model.
    """
    out=[];e=evidence();A=e['figure_7']['area_m2'];eps=.427
    for series,dp in [('calibration',.8),('holdout',.45)]:
        row=next(r for r in e['column_table_5'] if r['flow_mL_min']==43.8 and
            r['length_cm']==10 and r['inlet_mg_L']==150 and r['diameter_mm']==dp)
        t,y=observations(series);u=row['flow_mL_min']*1e-6/60/A
        Pe=row['derived_Pe'];Vmin=row['fitted_Vmin_cm']*.01
        for volume,velocity in [('total column area',u),('void area',u/eps)]:
            V=velocity*t
            for expression,prefactor in [('printed Pe/4',Pe/4),('square-root alternative',math.sqrt(Pe/4))]:
                arg=prefactor*(V-Vmin)/np.sqrt(V*Vmin)
                pred=.5*(1+np.array([math.erf(float(x)) for x in arg]))
                out.append(dict(series=series,volume_area=volume,prefactor=expression,
                    RMSE=float(np.sqrt(np.mean((pred-y)**2))),
                    role='published-fit reconstruction audit; not validation or calibration'))
    return out


def plot(base, cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.1),layout='constrained',sharey=True)
    for ax,series in zip(axes,['calibration','holdout']):
        t,y=observations(series)
        ax.plot(t/60,y,'o',mfc='none',ms=4,label='Published markers, digitized')
        for label,result in cases[series].items():
            ax.plot(np.array(result['time_s'])/60,np.array(result['outlet_concentration_kg_m3'])/.15,
                    label=label,lw=1.5)
        ax.axhline(.05,color='0.6',ls=':',lw=1)
        ax.set(xlabel='Time (min)',xlim=(0,300),ylim=(0,1.02),
               title='0.8 mm: fitted' if series=='calibration' else '0.45 mm: held out')
        ax.grid(alpha=.2);ax.legend(fontsize=8,loc='lower right')
    axes[0].set_ylabel('Outlet / inlet phenol concentration (1)')
    fig.suptitle('Date-stone carbon: testing transfer of a reduced uptake model',fontsize=12)
    fig.savefig(base/'comparison.png',dpi=180);plt.close(fig)


def main():
    started=time.perf_counter();e=evidence()
    run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=ROOT/'results/adsorption/alhamed2009/runs'/run_id;base.mkdir(parents=True)
    write(base/'evidence.json',e)
    initialization=time.perf_counter();module();initialization=time.perf_counter()-initialization
    files=[*Path(__file__).parent.glob('*.py'),PIXELS,
        ROOT/'literature/alhamed2009_protocol.md',ROOT/'julia/adsorption/AdsorptionColumn.jl',
        ROOT/'adsorption/bridge.py',ROOT/'src/pore2process_demo/bridge.py',
        ROOT/'uv.lock',ROOT/'julia/Project.toml',ROOT/'julia/Manifest.toml']
    sources={str(p.relative_to(ROOT)):sha(p) for p in files}
    for p in files:
        target=base/'source_snapshot'/p.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    write(base/'source_sha256.json',sources)
    print('Fitting 0.8 mm only...',flush=True)
    fitted=fit();write(base/'fit.json',fitted)
    rate=fitted['effective_ldf_s_inv'];D=fitted['dispersion_m2_s']
    print('Calibration', {k:v for k,v in fitted.items() if k!='history'},flush=True)
    factor=(.8/.45)**2
    specs={'calibration':{'Fitted reduced model':(rate,D)},
           'holdout':{'Frozen dispersion; rate ~ 1/d²':(rate*factor,D),
                      'Scenario: dispersion ~ d':(rate*factor,D*.45/.8)}}
    results={};metrics={};refinement={}
    for series,variants in specs.items():
        results[series]={};metrics[series]={}
        for label,(r,d) in variants.items():
            p=inputs(r,d);result=column(**p)
            if result['normalized_mass_residual']>1e-8:
                raise ValueError('Conservation screen failed.')
            results[series][label]=result;metrics[series][label]=compare(result,series)
            key=f'{series}-{len(results[series])}'
            write(base/f'{key}.json',dict(inputs=p,result=result,metrics=metrics[series][label]))
            if len(results[series])==1:
                changes={}
                for name,cells,dt in [('space',160,.25),('time',80,.125)]:
                    refined=column(**inputs(r,d,cells,dt));m=compare(refined,series)
                    changes[name]={k:abs(m[k]['predicted_s']/metrics[series][label][k]['predicted_s']-1)
                                   for k in ['t05','t50']}
                    changes[name]['normalized_mass_residual']=refined['normalized_mass_residual']
                    write(base/f'{series}-{name}.json',dict(inputs=inputs(r,d,cells,dt),result=refined,metrics=m))
                refinement[series]=changes
    # The model's calibration t05 is only about six seconds. Resolve that
    # short transient separately; the full-duration solver has a deliberate
    # 200000-step ceiling. No tolerance or solver guard is changed.
    early={};early_cases={}
    for name,cells,dt in [('base',80,.05),('space',160,.05),('time',80,.025)]:
        p=inputs(rate,D,cells,dt);p['t_end_s']=120.
        result=column(**p);m=compare(result,'calibration')
        early_cases[name]=m['t05']['predicted_s']
        early[name]=dict(t05_s=m['t05']['predicted_s'],normalized_mass_residual=result['normalized_mass_residual'])
        write(base/f'calibration-early-{name}.json',dict(inputs=p,result=result,metrics=m))
    early['space_relative']=abs(early_cases['space']/early_cases['base']-1)
    early['time_relative']=abs(early_cases['time']/early_cases['base']-1)
    plot(base,results)
    primary=next(iter(metrics['holdout'].values()))
    numerical=all(value<.02 for series,c in refinement.items() for r in c.values()
                  for key,value in r.items() if key=='t50' or key=='t05' and series=='holdout')
    numerical &= all(r['normalized_mass_residual']<=1e-8 for c in refinement.values() for r in c.values())
    numerical &= early['space_relative']<.02 and early['time_relative']<.02
    numerical &= all(early[k]['normalized_mass_residual']<=1e-8 for k in ['base','space','time'])
    calibration=next(iter(metrics['calibration'].values()))
    calibrated=calibration['RMSE_screen_passed'] and calibration['t50_screen_passed']
    accepted=numerical and calibrated and primary['RMSE_screen_passed'] and primary['t50_screen_passed'] and primary['t05_screen_passed'] is not False
    status=('passed_bounded_particle_size_comparison' if accepted else
            'failed_calibration_and_particle_size_transfer' if not calibrated else 'failed_particle_size_transfer')
    summary=dict(status=status,run_id=run_id,source_pdf_sha256=sha(PDF),source_sha256=sources,
        evidence=e,fit={k:v for k,v in fitted.items() if k!='history'},metrics=metrics,
        refinement=refinement,numerical_refinement_passed=numerical,
        early_calibration_t05_refinement=early,
        refinement_scope='Full-duration t50 for both curves and t05 for holdout; separate 120 s calculation for the six-second calibration t05.',
        calibration_passed=calibrated,physical_comparison_passed=accepted,
        source_equation_10_audit=source_equation_audit(),
        pressure_validation=False,new_OpenLB_runs=0,new_DWSIM_runs=0,
        versions=runtime_versions(),Julia_initialization_s=initialization,
        limitations=['Same-paper holdout, not independent experimental campaign.',
            'Digitized raster readings with unreported experimental uncertainty.',
            'Effective LDF fit does not separately identify film and pore diffusion.',
            'No native hydraulic specimen validation or full plant energy optimum.'],
        elapsed_s=time.perf_counter()-started)
    write(base/'summary.json',summary)
    with (base/'digitized_points.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['series','time_s','Cout_Cin','particle_diameter_m'])
        for series,p in points().items():
            writer.writerows((series,float(t),float(y),p['diameter_m']) for t,y in zip(p['time_s'],p['fraction']))
    summary['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*')
                                 if p.is_file() and p.name!='summary.json'}
    write(base/'summary.json',summary)
    write(ROOT/'results/adsorption/alhamed2009/latest.json',{'summary':str(base/'summary.json')})
    print(json.dumps(dict(status=status,numerical_refinement_passed=numerical,metrics=metrics,summary=str(base/'summary.json')),indent=2))


if __name__=='__main__':
    main()
