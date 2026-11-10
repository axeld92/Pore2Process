"""Run stock D3Q27/TRT; check fields independently before publishing evidence."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import math
from pathlib import Path
import subprocess
import tempfile
import shutil
from uuid import uuid4
import numpy as np
from coupling.run import sha
from hydraulic.run import write_json
from hydraulic.build import ROOT,verify_vendor
from pore2process_demo.bridge import runtime,runtime_versions
from process_filtration.run import read
from star_channel.geometry import polygon,write_stl,dimensions
from star_channel.native import analyze
from .build import main as build,BINARY


def inspect(directory):
    a=analyze(directory);r=a['raw']
    h=np.atleast_1d(np.genfromtxt(Path(directory)/'history.csv',delimiter=',',names=True))
    a['energy_convergence_relative_std']=float(np.std(h['mean_lattice_kinetic_energy'][-6:],ddof=1)/abs(np.mean(h['mean_lattice_kinetic_energy'][-6:])))
    a['mean_wall_traction_from_momentum_Pa']=abs(r['gradient_Pa_m'])*a['geometry']['open_area_m2']/a['geometry']['wetted_perimeter_m']
    a['viscous_observation_time']=r['steps']*r['dt_s']*(r['viscosity_Pa_s']/r['density_kg_m3'])/a['geometry']['hydraulic_diameter_m']**2
    return a


def solve(directory,n,gradient,*,tau=.8,axial=4,circle=False,warm=None,warm_scale=1.,max_steps=1800000):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    odd=.5+.1875/(tau-.5);dx=.0046/n;dt=(tau-.5)*dx*dx/(3e-6)
    Dh=.0046 if circle else dimensions(polygon())['hydraulic_diameter_m']
    minimum=math.ceil(.05*Dh**2/(1e-6*dt))
    args=[str(BINARY),'--resolution',str(n),'--gradient',str(gradient),'--tau',str(tau),
        '--odd-tau',str(odd),'--axial',str(axial),'--min-steps',str(minimum),'--max-steps',str(max_steps),'--snapshots']
    if circle:args+=['--circle']
    else:
        with tempfile.TemporaryDirectory(dir=ROOT/'.local/operating_map') as temp:
            expected=Path(temp)/'channel.stl';write_stl(expected,polygon(),-10*dx,(axial+10)*dx,.0023)
            target=directory/'channel.stl'
            if target.exists() and target.read_bytes()!=expected.read_bytes():raise ValueError('Resumed geometry changed.')
            if not target.exists():shutil.copy2(expected,target)
    if warm:args+=['--warm',str(Path(warm).resolve()),'--warm-scale',str(warm_scale)]
    request={'command':args,'binary_sha256':sha(BINARY),'warm_sha256':None if warm is None else sha(Path(warm)),
        'scope':'Native D3Q27/TRT, constant magic 3/16; initial estimate from measured native velocity only.'}
    request_path=directory/'request.json'
    if request_path.exists() and read(request_path)!=request:raise ValueError('Resume input or binary changed: '+str(directory))
    write_json(request_path,request)
    if (directory/'raw.json').exists() and read(directory/'raw.json').get('converged'):
        print('Reusing hash-matched native case '+directory.name,flush=True);return inspect(directory)
    if (directory/'history.csv').exists():
        archive=directory/'attempts'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        archive.mkdir(parents=True)
        for name in ('history.csv','stdout.log','raw.json','section.csv','planes.csv'):
            f=directory/name
            if f.exists():shutil.move(str(f),archive/name)
    write_json(directory/'execution_policy.json',{'wall_time_limit_s':7200,'stationarity_tolerance':1e-8,'minimum_viscous_observation_time':.05})
    with (directory/'stdout.log').open('w') as f:
        subprocess.run(args,cwd=directory,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=7200)
    return inspect(directory)


def checks(cases,reference,p):
    lim=p['native_acceptance'];out={}
    for name,a in cases.items():
        r=a['raw']
        out[name+'_stock_TRT27']=r['collision_model']=='TRT' and r['velocity_set']=='D3Q27'
        out[name+'_magic']=math.isclose((r['tau']-.5)*(r['odd_tau']-.5),.1875,rel_tol=1e-10)
        out[name+'_stationarity']=a['convergence_relative_std']<=1e-8*(1+1e-5) and a['energy_convergence_relative_std']<=1e-8*(1+1e-5)
        out[name+'_observation_time']=a['viscous_observation_time']>=.05
        out[name+'_Mach']=r['max_Mach']<=lim['Mach']
        out[name+'_density']=r['max_density_departure']<=lim['density_departure']
        out[name+'_mass']=a['plane_mass_relative_spread']<=lim['plane_mass_relative']
        out[name+'_transverse']=r['max_transverse_speed_m_s']/r['max_speed_m_s']<=lim['transverse_relative']
    C=lambda k:cases[k]['geometry_coefficient_m4']
    out['star_grid']=abs(C('grid64')/C('grid48')-1)<=lim['grid_relative']
    out['star_independent_Julia']=abs(C('grid64')/reference-1)<=lim['reference_relative']
    out['circle_shape']=cases['circle64']['raw']['shape']=='circle' and cases['circle_high']['raw']['shape']=='circle'
    out['circle_analytical']=abs(C('circle64')/(np.pi*(.0046/2)**4/8)-1)<.02
    out['circle_high_pressure_law']=abs(C('circle_high')/C('circle32')-1)<.02
    out['high_Re_grid']=abs(C('high1388_n64')/C('high1388_n48')-1)<=lim['grid_relative']
    out['high_Re_time_step']=abs(C('time')/C('high1388_n64')-1)<=lim['time_step_relative']
    out['high_Re_axial_domain']=abs(C('domain')/C('high1388_n64')-1)<=lim['domain_relative']
    out['high_Re_initial_estimate']=abs(C('seed')/C('high1388_n32')-1)<1e-5
    for name,a in cases.items():
        if name.startswith('high'):
            out[name+'_pressure_law']=abs(a['geometry_coefficient_m4']/C('grid'+str(a['raw']['resolution']))-1)<=lim['pressure_scaling_relative']
    if not all(out.values()):raise ValueError('Native verification failed: '+str({k:v for k,v in out.items() if not v}))
    return out


def execute(output,resume=None):
    p=read(ROOT/'operating_map/protocol.json');meta=build();jl=runtime()
    base=Path(resume) if resume else Path(output)/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8])
    base=base.resolve();base.mkdir(parents=True,exist_ok=True)
    fingerprint={str(f.relative_to(ROOT)):sha(f) for f in (ROOT/'operating_map/openlb').iterdir() if f.is_file()}
    state={'protocol':p,'source_sha256':fingerprint,'binary_sha256':sha(BINARY)}
    if (base/'summary.json').exists():raise ValueError('Accepted runs are immutable; create a new run.')
    if (base/'inputs.json').exists() and read(base/'inputs.json')!=state:raise ValueError('Cannot resume changed native configuration.')
    write_json(base/'inputs.json',state)
    previous=Path(read(ROOT/'results/star_channel/native/latest.json')['summary']);old=read(previous)
    cases={}
    for n in (32,48,64):
        print('Native baseline grid'+str(n),flush=True)
        cases['grid'+str(n)]=solve(base/('grid'+str(n)),n,1.,warm=previous.parent/f'grid{n}/section.csv')
    print('Native circular oracle',flush=True)
    cases['circle64']=solve(base/'circle64',64,1.,circle=True)
    cases['circle32']=solve(base/'circle32',32,1.,circle=True)
    if (base/'cold64/raw.json').exists():cases['cold64']=inspect(base/'cold64')
    jobs=[]
    for re in (1000,1388,1600,1850,2000):
        jobs.append((f'high{re}_n64',64,re,{}))
    jobs += [('high1388_n48',48,1388,{}),('time',64,1388,{'tau':.5005}),
             ('domain',64,1388,{'axial':8}),('seed',32,1388,{'warm_fraction':.98}),
             ('circle_high',32,1388,{'circle':True}),('high1388_n32',32,1388,{})]
    def run(job):
        name,n,re,kw=job;kw=dict(kw);circle=kw.pop('circle',False);fraction=kw.pop('warm_fraction',1.)
        baseline=cases[('circle' if circle else 'grid')+str(n)]
        gradient=re/baseline['Re_mean'];tau=kw.pop('tau',.5+n/64000)
        print('Native operating case '+name,flush=True)
        initial=base/(('circle' if circle else 'grid')+str(n))/'section.csv';scale=gradient*fraction
        restart=base/'restart_initial_fields'/name/'section.csv'
        if restart.exists():
            initial=restart;scale=1.
            parent=read(restart.parent/'source_request.json')
            if parent['binary_sha256']!=sha(BINARY):raise ValueError('Restart estimate came from a different native binary.')
        return name,solve(base/name,n,gradient,tau=tau,circle=circle,
            warm=initial,warm_scale=scale,**kw)
    # Start the longest independent circular case first.
    jobs.sort(key=lambda j:j[0]!='circle_high')
    with ThreadPoolExecutor(max_workers=4) as pool:
        for name,value in pool.map(run,jobs):cases[name]=value
    vertices=jl.seval('Matrix{Float64}')(polygon());ref=float(jl.duct_stokes_reference(vertices,256).conductance_geometry_m4)
    accepted=checks(cases,ref,p)
    if fingerprint!={str(f.relative_to(ROOT)):sha(f) for f in (ROOT/'operating_map/openlb').iterdir() if f.is_file()}:raise ValueError('Native source changed during verification.')
    result={'status':'passed',**state,'build':meta,'cases':cases,'accepted_case':'high1388_n64',
        'geometry_coefficient_m4':cases['high1388_n64']['geometry_coefficient_m4'],
        'verified_Re_range':[cases['high1000_n64']['Re_mean'],cases['high2000_n64']['Re_mean']],
        'geometry':dimensions(polygon()),'reference_C_m4':ref,'checks':accepted,'versions':runtime_versions(),
        'scope':'Fully developed laminar branch of a straight periodic nonpermeating star. No entry loss, turbulent-transition or permeating-wall CFD validation.',
        'original_vendor':verify_vendor()}
    result['artifacts_sha256']={str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file()}
    write_json(base/'summary.json',result);write_json(Path(output)/'latest.json',{'summary':str(base/'summary.json')})
    print('Native envelope accepted: '+str(base),flush=True);return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/operating_map/native')
    parser.add_argument('--resume',type=Path);a=parser.parse_args();execute(a.output,a.resume)
