"""Compare saved stock-native fields against the independent continuum solution."""
import json
from pathlib import Path

import numpy as np

from . import reference as ref

FIELD_NAMES = ('initial','middle','peak','final')
SHELL_M = 50e-6
RMS_LIMIT = .01
EXCURSION_LIMIT = .02


def load(directory):
    directory=Path(directory);r=json.loads((directory/'raw.json').read_text())
    shape=(r['n'],r['n'],r['nz']);h=r['dx_m']
    mask=np.fromfile(directory/'material.u8',dtype='u1').reshape(shape)
    u=np.fromfile(directory/'velocity.f64',dtype='<f8').reshape((*shape,3))
    f={k:np.fromfile(directory/f'{k}.f64',dtype='<f8').reshape(shape) for k in FIELD_NAMES}
    history=np.genfromtxt(directory/'inventory.csv',delimiter=',',names=True)
    q=[(np.arange(n)+.5)*h for n in shape]
    return r,mask,u,f,history,q


def profile(kind,q,time,speed,amplitude,mode=None):
    x,y,z=q
    if kind=='disk':return ref.disk_profile(x[:,None,None],y[None,:,None],time,speed,amplitude)+np.zeros_like(z)[None,None,:]
    if kind=='uniform':
        return (.5+amplitude*np.exp(-ref.D*ref.K**2*time)*np.sin(ref.K*(x-speed*time)))[:,None,None]+np.zeros((len(y),len(z)))
    f=1. if mode is None else ref.sample(mode,y,z)
    return .5+amplitude*np.imag(np.exp(1j*ref.K*x[:,None,None])*f)+np.zeros((len(y),len(z)))


def analyze(directory,modes=None):
    directory=Path(directory);r,mask,u,f,h,q=load(directory)
    request=json.loads((directory/'request.json').read_text())['spec'];kind=r['kind'];speed=r['speed_m_s'];amp=r['amplitude_kg_m3']
    steps=r['steps'];dt=r['dt_s'];dx=r['dx_m'];disk=kind=='disk';fluid=mask==1
    if (kind,r['n'],r['nz'],steps,dt,speed,amp)!=(request['kind'],request['n'],request['nz'],request['steps'],request['dt'],request['speed'],request['amplitude']):
        raise ValueError('Native parameters differ from the recorded request')
    if r['side_m']!=ref.L or r['diffusivity_m2_s']!=ref.D or r['radius_m']!=ref.R:
        raise ValueError('Physical scales differ from the benchmark')
    if (r['ax'],r['beta'])!=(ref.AX,ref.B) or 'AdvectionDiffusionTRTdynamics' not in r['collision']:
        raise ValueError('Bulk flow/collision contract differs')
    for v in (*f.values(),u,*(h[k] for k in h.dtype.names)):
        if not np.isfinite(v).all():raise ValueError('Non-finite saved field/history')
    np.testing.assert_array_equal(h['step'],np.arange(steps+1))
    np.testing.assert_allclose(h['time_s'],np.arange(steps+1)*dt,rtol=1e-14,atol=0)
    x,y,z=np.meshgrid(*q,indexing='ij',sparse=False)
    radial=np.hypot(x-ref.L/2,y-ref.L/2)
    expected_mask=np.where(radial<=ref.R,1,2) if disk else np.ones_like(mask)
    np.testing.assert_array_equal(mask,expected_mask)
    expected_velocity=ref.velocity(kind,x,y,z,speed);expected_velocity[~fluid]=0
    velocity_error=float(np.max(abs(u-expected_velocity)))
    initial_error=float(np.max(abs(f['initial']-profile(kind,q,0.,speed,amp))))
    if velocity_error>1e-14 or initial_error>1e-12:raise ValueError('Incorrect saved velocity or initialized field')
    restored_odd=r['magic']/(r['tau_even']-.5)
    recovered_D=restored_odd*dx*dx/(3*dt)
    if abs(recovered_D/ref.D-1)>1e-9:raise ValueError('TRT odd relaxation changes molecular diffusion')
    max_lattice_speed=float(np.max(np.linalg.norm(u,axis=-1)))*dt/dx
    if max_lattice_speed>.05:raise ValueError('Lattice speed outside declared range')
    # Verify the reported extrema against the full saved snapshots.
    for key,step in [('initial',0),('middle',steps//2),('final',steps),('peak',r['peak_step'])]:
        np.testing.assert_allclose([np.min(f[key][fluid]),np.max(f[key][fluid])],
            [h['min_concentration_kg_m3'][step],h['max_concentration_kg_m3'][step]],rtol=0,atol=2e-14)
    lower=float(np.min(h['min_concentration_kg_m3']));upper=float(np.max(h['max_concentration_kg_m3']))
    excursion=max(.5-amp-lower,upper-(.5+amp),0.)
    if abs(excursion-r['maximum_excursion_kg_m3'])>1e-14:raise ValueError('Peak and history disagree')
    fluid_mass=h['fluid_inventory_kg'];kin=h['kinetic_inventory_kg'];volume=int(np.count_nonzero(fluid))*dx**3
    if abs(fluid_mass[0]/(.5*volume)-1)>1e-11:raise ValueError('Initial liquid inventory differs from voxel volume')
    mass_drift=float(np.max(abs(kin/kin[0]-1)));liquid_drift=float(np.max(abs(fluid_mass/fluid_mass[0]-1)))
    if mass_drift>1e-8:raise ValueError('Closed numerical inventory fails its declared screen')
    adjacent=np.zeros_like(fluid);links=0;expected_wall_mass=0.
    if disk:
        stencil=np.genfromtxt(directory/'stencil.csv',delimiter=',',names=True)
        for row in stencil:
            direction=tuple(int(row[k]) for k in ('cx','cy','cz'))
            incoming=np.roll(fluid,direction,axis=(0,1,2)) & ~fluid
            links+=int(np.count_nonzero(incoming));expected_wall_mass+=float(row['weight']*np.sum(f['initial'][incoming]))*dx**3
            adjacent|=np.roll(incoming,tuple(-i for i in direction),axis=(0,1,2))
        if links!=r['active_wall_links'] or abs(expected_wall_mass/h['active_wall_population_kg'][0]-1)>1e-11:
            raise ValueError('Active-wall link count or initial storage disagrees')
        if max(float(np.max(abs(v-v[:,:,:1]))) for v in f.values())>1e-12:
            raise ValueError('Axial invariance failed in the disk')
    elif r['active_wall_links']!=0 or np.any(h['active_wall_population_kg']!=0):raise ValueError('Wall in the periodic cube')
    stats=dict(raw=r,range_kg_m3=[lower,upper],maximum_excursion_kg_m3=excursion,
        closed_inventory_relative_drift=mass_drift,physical_liquid_relative_drift=liquid_drift,
        initial_fluid_volume_m3=volume,velocity_maximum_error_m_s=velocity_error,
        initial_maximum_error_kg_m3=initial_error,max_lattice_speed=max_lattice_speed,
        reconstructed_diffusivity_m2_s=recovered_D,peak_time_s=r['peak_step']*dt,
        checks={'engineering_bounds':excursion<=EXCURSION_LIMIT,'physical_bounds':excursion<=1e-8})
    if disk:
        stats.update(analytic_fluid_volume_m3=np.pi*ref.R**2*r['nz']*dx,
                     initial_active_wall_inventory_kg=expected_wall_mass,active_wall_links=links,
                     wall_shell_m=SHELL_M,adjacent_fraction=float(np.mean(adjacent[fluid])))
    if amp==0:
        stats['checks']['constant']=excursion<=1e-9
        return stats
    stats['fields']={}
    for i,key in enumerate(('middle','final')):
        time=(steps//2 if i==0 else steps)*dt
        if kind=='stretch' and (modes is None or abs(time-ref.TIMES[i])>1e-14):raise ValueError('Wrong stretching reference time')
        exact=profile(kind,q,time,speed,amp,mode=None if modes is None else modes[i]);err=(f[key]-exact)/amp
        rms=lambda subset:float(np.sqrt(np.mean(err[subset]**2)))
        stats['fields'][key]=dict(time_s=time,RMS_over_amplitude=rms(fluid),max_over_amplitude=float(np.max(abs(err[fluid]))))
        stats['checks'][f'{key}_profile']=rms(fluid)<=RMS_LIMIT
        if disk:
            shell=fluid&(radial>=ref.R-SHELL_M)
            stats['fields'][key].update(wall_shell_RMS_over_amplitude=rms(shell),
                                        wall_link_RMS_over_amplitude=rms(adjacent),
                                        interior_RMS_over_amplitude=rms(fluid&~shell))
            stats['checks'][f'{key}_wall_shell_profile']=rms(shell)<=RMS_LIMIT
            cos=(ref.disk_profile(x,y,0.,speed,1.)-.5)[fluid]
            theta=np.arctan2(y-ref.L/2,x-ref.L/2)
            radial_mode=ref.bessel(2,ref.ROOT_J2_PRIME*radial/ref.R)/ref.bessel(2,ref.ROOT_J2_PRIME)
            sin=(radial_mode*np.sin(2*theta))[fluid]
            coefficients=np.linalg.lstsq(np.column_stack([cos,sin]),(f[key][fluid]-.5)/amp,rcond=None)[0]
            expected_decay=np.exp(-ref.D*(ref.ROOT_J2_PRIME/ref.R)**2*time)
            phase=np.angle(complex(*coefficients)*np.exp(-2j*speed/ref.R*time),deg=True)
            stats['fields'][key].update(mode_contrast=float(np.linalg.norm(coefficients)),
                analytical_mode_contrast=float(expected_decay),pattern_angle_error_deg=float(phase/2))
    return stats


def compare_depth(thin,full):
    a,am,au,af,ah,_=load(thin);b,bm,bu,bf,bh,_=load(full)
    difference={k:float(np.max(abs(af[k][:,:,:1]-bf[k]))) for k in af}
    factor=b['nz']/a['nz'];mass={}
    np.testing.assert_array_equal(bm,np.broadcast_to(am[:,:,:1],bm.shape))
    np.testing.assert_array_equal(bu,np.broadcast_to(au[:,:,:1,:],bu.shape))
    for k in ('min_concentration_kg_m3','max_concentration_kg_m3'):
        np.testing.assert_allclose(ah[k],bh[k],rtol=0,atol=1e-12)
    for k in ('fluid_inventory_kg','active_wall_population_kg','kinetic_inventory_kg'):
        mass[k]=float(np.max(abs(bh[k]/(ah[k]*factor)-1)))
    if max(difference.values())>1e-12 or max(mass.values())>1e-11:
        raise ValueError('Thin disk does not reproduce the full-depth result')
    return dict(field_max_abs_differences_kg_m3=difference,mass_scaling_relative_errors=mass,depth_ratio=factor,passed=True)
