"""Execute and inspect actual OpenLB fields; no Python hydraulic solver."""
import json
from pathlib import Path
import subprocess
import numpy as np

from hydraulic.run import write_json
from .build import BINARY
from .geometry import polygon, dimensions, write_stl

def solve(directory, resolution, *, gradient=1., axial=4, tau=.8, mu=.001, rho=1000.,
          circle=False, warm=None, warm_scale=1., max_steps=120000,binary=BINARY,odd_tau=1.6,snapshots=False):
    directory = Path(directory); directory.mkdir(parents=True)
    d = .0046/resolution
    if not circle:
        write_stl(directory/'channel.stl',polygon(),-10*d,(axial+10)*d,.0023)
    command = [str(binary),'--resolution',str(resolution),'--gradient',str(gradient),
               '--axial',str(axial),'--tau',str(tau),'--mu',str(mu),'--rho',str(rho),
               '--odd-tau',str(odd_tau),
               '--max-steps',str(max_steps)]
    if circle: command.append('--circle')
    if snapshots:command.append('--snapshots')
    if warm:
        command += ['--warm',str(Path(warm).resolve()),'--warm-scale',str(warm_scale)]
    write_json(directory/'request.json',{'command':command,'warm_start_native_field':None if warm is None else str(warm)})
    with (directory/'stdout.log').open('w') as log:
        subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=2400)
    return analyze(directory)

def analyze(directory):
    directory=Path(directory); r=json.loads((directory/'raw.json').read_text())
    section=np.atleast_1d(np.genfromtxt(directory/'section.csv',names=True,delimiter=','))
    planes=np.atleast_1d(np.genfromtxt(directory/'planes.csv',names=True,delimiter=','))
    history=np.atleast_1d(np.genfromtxt(directory/'history.csv',names=True,delimiter=','))
    for array in (section,planes,history):
        if not all(np.isfinite(array[n]).all() for n in array.dtype.names):
            raise ValueError('Non-finite native field.')
    if r['solver']!='OpenLB' or r['converged'] is not True or len(section)!=r['section_nodes']:
        raise ValueError('Changed native field contract.')
    if not np.allclose(planes['Q_m3_s'],r['Q_m3_s'],rtol=1e-11,atol=0):
        raise ValueError('Periodic-section flow does not match native JSON.')
    for name,coord in [('j','y_m'),('k','z_m')]:
        # The mesh starts half a cell before y=z=0; indices here are OpenLB
        # local block indices, explicitly distinct from Python export node IDs.
        if not np.allclose(section[coord],(section[name]-.5)*r['dx_m'],rtol=0,atol=1e-14):
            raise ValueError('Native physical-coordinate/index convention changed.')
    std=float(np.std(history['Q_m3_s'][-6:],ddof=1)/abs(np.mean(history['Q_m3_s'][-6:])))
    if std>r['tolerance']*(1+1e-5): raise ValueError('Exported convergence history fails native criterion.')
    C=r['viscosity_Pa_s']*r['Q_m3_s']/r['gradient_Pa_m']
    dims=dimensions(polygon()) if r['shape']=='six_point_star' else {
        'open_area_m2':np.pi*(r['diameter_m']/2)**2,'wetted_perimeter_m':np.pi*r['diameter_m'],
        'hydraulic_diameter_m':r['diameter_m']}
    area_voxel=len(section)*r['dx_m']**2
    return {'raw':r,'geometry':dims,'geometry_coefficient_m4':C,
        'mean_velocity_m_s':r['Q_m3_s']/dims['open_area_m2'],
        'Re_mean':abs(r['Q_m3_s'])/dims['open_area_m2']*dims['hydraulic_diameter_m']*r['density_kg_m3']/r['viscosity_Pa_s'],
        'plane_mass_relative_spread':float(np.ptp(planes['mass_kg_s'])/abs(np.mean(planes['mass_kg_s']))),
        'convergence_relative_std':std,'voxel_area_m2':area_voxel,
        'voxel_area_relative_error':area_voxel/dims['open_area_m2']-1}

def verify_suite(cases,references,p):
    lim=p['acceptance']; checks={}
    for name,c in cases.items():
        checks[name+'_mass']=c['plane_mass_relative_spread']<=lim['plane_mass_relative']
        checks[name+'_Mach']=c['raw']['max_Mach']<=lim['Mach']
        checks[name+'_collision']=c['raw'].get('collision_model')==p['native']['collision_model']
        checks[name+'_odd_tau']=bool(np.isclose(c['raw'].get('odd_tau',0),p['native']['odd_tau']))
    C=lambda name:cases[name]['geometry_coefficient_m4']
    finest=f"grid{p['native'].get('accepted_resolution',48)}"
    previous=f"grid{p['native']['resolutions'][-2]}"
    checks['grid']=abs(C(finest)/C(previous)-1)<=lim['grid_relative']
    checks['independent_reference']=abs(C(finest)/references[-1]['conductance_geometry_m4']-1)<=lim['independent_reference_relative']
    checks['reference_grid']=abs(references[-1]['conductance_geometry_m4']/references[-2]['conductance_geometry_m4']-1)<.01
    exact=np.pi*(.0046/2)**4/8
    checks['circle_analytical']=abs(C('circle')/exact-1)<.02
    for name in ('half','double','reverse'):
        checks[name+'_linear']=abs(C(name)/C('grid32')-1)<=lim['scaling_relative']
    checks['periodic_domain']=abs(C('domain')/C('grid32')-1)<=lim['domain_relative']
    checks['relaxation_time']=abs(C('tau')/C('grid32')-1)<=lim['tau_relative']
    if 'op20_n64' in cases:
        for name,c in cases.items():
            if name.startswith('op'):
                reference=f"grid{c['raw']['resolution']}"
                checks[name+'_pressure_law']=abs(C(name)/C(reference)-1)<=lim['scaling_relative']
                checks[name+'_transverse']=c['raw']['max_transverse_speed_m_s']/c['raw']['max_speed_m_s']<.01
        checks['operating_grid']=abs(C('op20_n64')/C('op20_n48')-1)<=lim['grid_relative']
        checks['operating_domain']=abs(C('op20_domain')/C('op20_n64')-1)<=lim['domain_relative']
        checks['operating_time_step']=abs(C('op20_time')/C('op20_n64')-1)<=lim['tau_relative']
        checks['operating_initial_condition']=abs(C('op20_cold')/C('op20_n48')-1)<1e-7
    if not all(checks.values()): raise ValueError(f'Native verification failed: {checks}')
    return checks
