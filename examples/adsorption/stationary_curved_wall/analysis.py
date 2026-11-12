"""Saved native data checks and comparison against the radial reference."""
import json
from pathlib import Path

import numpy as np

from examples.adsorption.transport_verification.analysis import load,SHELL_M,RMS_LIMIT,EXCURSION_LIMIT
from . import reference as ref


def geometry(directory,fluid,u,initial,h):
    stencil=np.genfromtxt(directory/'stencil.csv',delimiter=',',names=True)
    adjacent=np.zeros_like(fluid);count=0;wall=0.;normal=[]
    for row in stencil:
        p=tuple(int(row[k]) for k in ('cx','cy','cz'))
        incoming=np.roll(fluid,p,axis=(0,1,2)) & ~fluid
        origin=np.roll(incoming,tuple(-i for i in p),axis=(0,1,2))
        count+=int(np.count_nonzero(incoming));adjacent|=origin
        wall+=float(row['weight']*np.sum(initial[incoming]))*h**3
        if np.any(origin):normal.append(u[origin]@np.asarray(p)/np.linalg.norm(p))
    link_speed=np.concatenate(normal)
    return adjacent,count,wall,dict(maximum_adjacent_speed_m_s=float(np.max(np.linalg.norm(u[adjacent],axis=-1))),
        RMS_link_normal_speed_m_s=float(np.sqrt(np.mean(link_speed**2))),
        maximum_absolute_link_normal_speed_m_s=float(np.max(abs(link_speed))))


def rms(x):return float(np.sqrt(np.mean(x*x)))


def analyze(directory,radial):
    directory=Path(directory);r,mask,u,fields,h,q=load(directory)
    req=json.loads((directory/'request.json').read_text())['spec'];a=req['amplitude']
    n=r['n'];nz=r['nz'];dx=r['dx_m'];dt=r['dt_s'];steps=r['steps']
    if (r['kind'],n,nz,steps,dt,r['amplitude_kg_m3'],r['speed_m_s'])!=('stationary_disk',req['n'],req['nz'],ref.STEPS,ref.DT,a,ref.U):
        raise ValueError('Native case differs from request')
    if (r['side_m'],r['radius_m'],r['diffusivity_m2_s'])!=(ref.L,ref.R,ref.D) or 'AdvectionDiffusionTRTdynamics' not in r['collision']:
        raise ValueError('Changed physical scales or scalar dynamics')
    for v in (u,*fields.values(),*(h[k] for k in h.dtype.names)):
        if not np.isfinite(v).all():raise ValueError('Non-finite native data')
    np.testing.assert_array_equal(h['step'],np.arange(steps+1))
    np.testing.assert_allclose(h['time_s'],np.arange(steps+1)*dt,rtol=1e-14,atol=0)
    x,y=np.meshgrid(q[0],q[1],indexing='ij');rho=np.hypot(x-ref.L/2,y-ref.L/2)
    x=x[:,:,None];y=y[:,:,None];rho=rho[:,:,None];fluid=mask==1
    np.testing.assert_array_equal(mask,np.broadcast_to(np.where(rho<=ref.R,1,2),mask.shape))
    ux,uy=ref.velocity(x,y)
    velocity_error=max(float(np.max(abs(u[...,0]-np.where(fluid,ux,0.)))),
                       float(np.max(abs(u[...,1]-np.where(fluid,uy,0.)))),float(np.max(abs(u[...,2]))))
    initial_error=float(np.max(abs(fields['initial']-ref.profile(x,y,0.,a))))
    recovered=r['magic']/(r['tau_even']-.5)*dx**2/(3*dt)
    max_lattice_speed=float(np.max(np.linalg.norm(u,axis=-1)))*dt/dx
    if max(velocity_error/1e-14,initial_error/1e-12,abs(recovered/ref.D-1)/1e-9,max_lattice_speed/.05)>1:
        raise ValueError('Velocity, initialization or TRT diffusivity mismatch')
    theta=np.linspace(0.,2*np.pi,129)
    wall_speed=float(np.max(np.hypot(*ref.velocity(ref.L/2+ref.R*np.cos(theta),ref.L/2+ref.R*np.sin(theta)))))
    if wall_speed>1e-15:raise ValueError('Continuum velocity does not vanish at wall')
    z_error=max(float(np.max(abs(f-f[:,:,:1]))) for f in fields.values())
    if z_error>1e-12:raise ValueError('Disk field is not axially invariant')
    for key,index in [('initial',0),('middle',steps//2),('peak',r['peak_step']),('final',steps)]:
        f=fields[key][fluid]
        np.testing.assert_allclose([np.min(f),np.max(f)],[h['min_concentration_kg_m3'][index],h['max_concentration_kg_m3'][index]],rtol=0,atol=2e-14)
    low=float(np.min(h['min_concentration_kg_m3']));high=float(np.max(h['max_concentration_kg_m3']))
    excursion=max(.5-a-low,high-.5-a,0.)
    if abs(excursion-r['maximum_excursion_kg_m3'])>1e-14:raise ValueError('Peak excursion does not match history')
    fluid_inventory=h['fluid_inventory_kg'];kinetic=h['kinetic_inventory_kg'];volume=int(np.count_nonzero(fluid))*dx**3
    adjacent,links,wall,velocity_links=geometry(directory,fluid,u,fields['initial'],dx)
    if links!=r['active_wall_links'] or abs(wall/h['active_wall_population_kg'][0]-1)>1e-11:
        raise ValueError('Active wall storage check failed')
    mass=float(np.max(abs(kinetic/kinetic[0]-1)))
    if mass>1e-8 or abs(fluid_inventory[0]/(.5*volume)-1)>1e-11:raise ValueError('Numerical or initial liquid inventory failed')
    peak=fields['peak'];ijk=np.unravel_index(np.argmax(np.where(fluid,abs(peak-.5),-np.inf)),peak.shape)
    result=dict(raw=r,range_kg_m3=[low,high],maximum_excursion_kg_m3=excursion,peak_time_s=r['peak_step']*dt,
        peak_point=dict(ijk=[int(i) for i in ijk],concentration_kg_m3=float(peak[ijk]),
                        smooth_wall_gap_m=float(ref.R-rho[ijk[0],ijk[1],0]),touches_wall_link=bool(adjacent[ijk])),
        closed_inventory_relative_drift=mass,liquid_inventory_relative_drift=float(np.max(abs(fluid_inventory/fluid_inventory[0]-1))),
        initial_liquid_volume_m3=volume,continuum_liquid_volume_m3=float(np.pi*ref.R**2*nz*dx),
        initial_maximum_error_kg_m3=initial_error,velocity_maximum_error_m_s=velocity_error,
        reconstructed_diffusivity_m2_s=recovered,continuum_maximum_wall_speed_m_s=wall_speed,
        max_lattice_speed=max_lattice_speed,axial_invariance_maximum_error=z_error,**velocity_links,
        checks=dict(engineering_bounds=excursion<=EXCURSION_LIMIT,physical_bounds=excursion<=1e-8))
    if not a:
        result['checks']['constant']=excursion<=1e-9
        return result
    s,modes=radial;shell=fluid&(rho>=ref.R-SHELL_M)
    result['fields']={}
    for key,t,f in zip(('middle','final'),ref.TIMES,modes):
        exact=np.broadcast_to(ref.profile(x,y,t,a,(s,f)),mask.shape)
        error=(fields[key]-exact)/a
        result['fields'][key]=dict(time_s=t,RMS_over_amplitude=rms(error[fluid]),
            shell_RMS_over_amplitude=rms(error[shell]),link_RMS_over_amplitude=rms(error[adjacent]),
            interior_RMS_over_amplitude=rms(error[fluid&~shell]),max_over_amplitude=float(np.max(abs(error[fluid]))),
            liquid_sd_kg_m3=float(np.std(fields[key][fluid])),reference_sd_at_cells_kg_m3=float(np.std(exact[fluid])))
        result['checks'][key+'_profile']=rms(error[fluid])<=RMS_LIMIT
        result['checks'][key+'_shell']=rms(error[shell])<=RMS_LIMIT
    return result


def compare_background(wave_dir,constant_dir,radial):
    rw,mask,u,w,_,q=load(wave_dir);rc,cm,cu,c,_,_=load(constant_dir)
    for k in ('n','nz','dt_s','steps','speed_m_s','diffusivity_m2_s'):
        if rw[k]!=rc[k]:raise ValueError('Background controls must be matched')
    np.testing.assert_array_equal(mask,cm);np.testing.assert_array_equal(u,cu)
    x=q[0][:,None,None];y=q[1][None,:,None];rho=np.hypot(x-ref.L/2,y-ref.L/2)
    fluid=mask==1;shell=fluid&(rho>=ref.R-SHELL_M);result={};s,modes=radial
    for key,t,f in zip(('middle','final'),ref.TIMES,modes):
        original=(w[key]-ref.profile(x,y,t,ref.A,(s,f)))/ref.A
        background=(c[key]-.5)/ref.A
        result[key]={name:dict(original_RMS_over_amplitude=rms(original[region]),
                        background_RMS_over_amplitude=rms(background[region]),
                        diagnostic_wave_only_RMS_over_amplitude=rms((original-background)[region]))
                     for name,region in [('liquid',fluid),('outer_shell',shell)]}
    return dict(fields=result,production_correction_applied=False)
