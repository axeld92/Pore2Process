"""Deterministic bounded searches; all forward physics runs in Julia."""
import math

import numpy as np

from .bridge import batch, bed, call
from .cases import batch_inputs, column_inputs
from .data import observations, evidence


def golden(function, lo, hi, iterations=36):
    ratio=(math.sqrt(5)-1)/2
    a=hi-ratio*(hi-lo);b=lo+ratio*(hi-lo)
    fa=function(a);fb=function(b)
    for _ in range(iterations):
        if fa<fb:
            hi=b;b=a;fb=fa;a=hi-ratio*(hi-lo);fa=function(a)
        else:
            lo=a;a=b;fa=fb;b=lo+ratio*(hi-lo);fb=function(b)
    candidates=[(a,fa),(b,fb),(lo,function(lo)),(hi,function(hi))]
    return min(candidates,key=lambda x:x[1])


def endpoint_audit():
    rows=evidence()['batch_table_3']
    records=[]
    for row in rows:
        p=dict(initial_concentration_kg_m3=row['inlet_mg_L']*.001,
            fluid_volume_m3=row['solution_volume_L']*.001,adsorbent_mass_kg=row['dose_g']*.001,
            capacity_kg_kg=.1115,affinity_m3_kg=139.,linear_partition_m3_kg=None)
        result=call('bath_equilibrium',**p)
        records.append(dict(experiment=row['experiment'],inputs=p,result=result,
            reported_fitted_loading_kg_kg=row['fitted_qe_mg_g']*.001,
            relative_error=abs(result['loading_kg_kg']/(row['fitted_qe_mg_g']*.001)-1),
            role='training endpoint' if row['experiment'] in [4,5,6] else 'withheld endpoint condition'))
    def objective(factor):
        errors=[]
        for r in records:
            if r['experiment'] not in [4,5,6]:
                continue
            p={**r['inputs'],'capacity_kg_kg':.1115*factor}
            q=call('bath_equilibrium',**p)['loading_kg_kg']
            errors.append((q/r['reported_fitted_loading_kg_kg']-1)**2)
        return float(np.mean(errors))
    factor,error=golden(objective,.5,1.5)
    return dict(source_endpoints=records,source_max_relative_error=max(r['relative_error'] for r in records),
        source_endpoint_screen_passed=all(r['relative_error']<=.05 for r in records),
        diagnostic_capacity_multiplier=factor,diagnostic_training_relative_RMSE=math.sqrt(error),
        diagnostic_factor_applied_to_forward_model=False,primary_capacity_multiplier=1.,primary_affinity_multiplier=1.)


def fit_batch():
    selected={k:v for k,v in observations().items() if v['role']=='calibration'}
    history=[];cache={};bounds=np.array([[-18.,-10.],[-9.,-3.]])

    def objective(vector):
        vector=np.clip(vector,bounds[0],bounds[1]);key=tuple(np.round(vector,9))
        if key in cache:
            return cache[key]
        D,kf=10.**vector
        errors=[];residuals=[];reason=None
        try:
            for data in selected.values():
                p=batch_inputs(data['diameter_m'],D,kf,radial_cells=12,dt_s=120.)
                result=batch(**p)
                if result['normalized_mass_residual']>1e-8 or result['normalized_interface_flux_mismatch']>1e-8:
                    raise ValueError('Candidate failed numerical conservation/interface screen.')
                prediction=np.interp(data['time_s'],result['time_s'],result['mean_loading_kg_kg'])
                errors.extend((prediction-data['loading_kg_kg']).tolist())
                residuals.append(result['normalized_mass_residual'])
            score=float(np.mean(np.square(errors)))
        except Exception as error:
            reason=str(error).split('Stacktrace:')[0].strip()
            if not any(word in reason.lower() for word in ['converg','line search','negative','numerical conservation']):
                raise
            score=1e3
        history.append(dict(log10_Ds=float(vector[0]),log10_batch_kf_ref=float(vector[1]),
            MSE=None if reason else score,numerical_rejection=reason,
            max_normalized_mass_residual=max(residuals,default=None)))
        cache[key]=score
        return score

    seeds=[np.array([d,k]) for d in [-18,-16,-14,-12,-11,-10,-9] for k in [-10,-8,-6,-4,-3]]
    def search(start):
        simplex=[start.astype(float),np.clip(start+np.array([.5,0.]),bounds[0],bounds[1]),
                 np.clip(start+np.array([0.,.5]),bounds[0],bounds[1])]
        if np.array_equal(simplex[1],simplex[0]):simplex[1]=start-np.array([.5,0.])
        if np.array_equal(simplex[2],simplex[0]):simplex[2]=start-np.array([0.,.5])
        converged=False
        for iteration in range(240):
            simplex=sorted(simplex,key=objective);values=[objective(x) for x in simplex]
            if max(np.linalg.norm(x-simplex[0]) for x in simplex)<.01:
                converged=True;break
            centre=(simplex[0]+simplex[1])/2
            reflected=np.clip(2*centre-simplex[2],bounds[0],bounds[1]);fr=objective(reflected)
            if fr<values[0]:
                expanded=np.clip(centre+2*(reflected-centre),bounds[0],bounds[1])
                simplex[2]=expanded if objective(expanded)<fr else reflected
            elif fr<values[1]:
                simplex[2]=reflected
            else:
                contracted=np.clip(centre+.5*((reflected if fr<values[2] else simplex[2])-centre),bounds[0],bounds[1])
                if objective(contracted)<min(fr,values[2]):
                    simplex[2]=contracted
                else:
                    simplex=[simplex[0],(simplex[0]+simplex[1])/2,(simplex[0]+simplex[2])/2]
        best=min(simplex,key=objective)
        if not converged or objective(best)>=1e3:
            raise RuntimeError('Batch transport search did not converge to a numerically admitted candidate.')
        return best

    best=search(min(seeds,key=objective))
    profiles=[]
    for fixed_k in [-10,-9,-8,-7,-6,-5,-4,-3]:
        x,error=golden(lambda d:objective(np.array([d,float(fixed_k)])),-18.,-9.,iterations=24)
        profiles.append(dict(log10_batch_kf_ref=fixed_k,best_Ds_m2_s=10**x,RMSE_mg_g=math.sqrt(error)*1000))
    profile_best=min(profiles,key=lambda x:x['RMSE_mg_g'])
    restart=np.array([math.log10(profile_best['best_Ds_m2_s']),profile_best['log10_batch_kf_ref']])
    alternative=search(restart)
    best=min([best,alternative],key=objective);score=objective(best)
    return dict(diffusivity_m2_s=float(10.**best[0]),batch_film_ref_m_s=float(10.**best[1]),
        reference_diameter_m=.0008,batch_film_diameter_exponent=-1.,
        fitted_diameter_ids=list(selected),log10_parameters=best.tolist(),
        bounds_log10=bounds.tolist(),converged=True,multistart='coarse grid and best independent kf profile',
        search_RMSE_mg_g=math.sqrt(score)*1000,search_radial_cells=12,search_dt_s=120.,
        at_bound=bool(np.any(np.minimum(best-bounds[0],bounds[1]-best)<.011)),
        kf_profiles=profiles,profile_role='deterministic trade-off scan; not statistical confidence bounds',
        history=history,evaluations=len(history))


def fit_column_film(diffusivity):
    from adsorption.alhamed.study import observations as column_observations
    t,y=column_observations('calibration');history=[]
    def objective(logfactor):
        p,_=column_inputs(.0008,diffusivity,multiplier=10**logfactor,
                           axial_cells=40,radial_cells=12,dt_s=5.)
        r=bed(**p)
        pred=np.interp(t,r['time_s'],np.array(r['outlet_concentration_kg_m3'])/.15)
        score=float(np.mean((pred-y)**2))
        history.append(dict(log10_multiplier=logfactor,MSE=score,
            normalized_mass_residual=r['normalized_mass_residual']))
        return score
    # Coarse global scan before a local bounded one-dimensional refinement.
    grid=np.linspace(-4,2,25);scores=[objective(float(x)) for x in grid]
    i=int(np.argmin(scores));lo=float(grid[max(0,i-1)]);hi=float(grid[min(len(grid)-1,i+1)])
    x,error=golden(objective,lo,hi,iterations=28)
    return dict(multiplier=10**x,RMSE=math.sqrt(error),fitted_curve='Figure 7, 0.8 mm only',
        frozen_diffusivity_m2_s=diffusivity,equilibrium_multiplier=1.,inventory_multiplier=1.,
        bounds_log10=[-4.,2.],at_bound=min(x+4,2-x)<.01,history=history)
