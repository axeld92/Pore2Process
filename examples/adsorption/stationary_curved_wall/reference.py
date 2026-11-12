"""Python interface and independent convergence checks for the Julia oracle."""
import json
from pathlib import Path
import subprocess

import numpy as np

from pore2process_demo.io import sha,write_json
from pore2process_demo.paths import ROOT
from examples.adsorption.transport_verification.reference import A,D,DT,L,R,STEPS,TIMES,U,disk_profile

HERE=Path(__file__).parent
OMEGA0=3*np.sqrt(3)*U/(2*R)


def velocity(x,y,speed=U):
    x=np.asarray(x)-L/2;y=np.asarray(y)-L/2
    q=(x*x+y*y)/R**2
    angular=3*np.sqrt(3)*speed/(2*R)*np.maximum(0.,1-q)
    return -angular*y,angular*x


def runtime():
    data=json.loads((ROOT/'.local/julia_runtime.json').read_text())
    if Path(data['project']).resolve()!=(ROOT/'julia').resolve():
        raise ValueError('Select the project Julia environment with bootstrap_julia.py')
    return [data['executable'],'--startup-file=no',f"--project={data['project']}"]


def read_mode(path):
    a=np.genfromtxt(path,delimiter=',',names=True)
    f=np.array([a['middle_real']+1j*a['middle_imag'],a['final_real']+1j*a['final_imag']])
    return a['r_over_R'],f


def sample(s,f,radius_fraction):
    return np.interp(radius_fraction,np.r_[0.,s,1.],np.r_[0j,f,f[-1]])


def profile(x,y,t,amplitude=A,mode=None):
    if mode is None:
        if t!=0 and amplitude:raise ValueError('A nonzero time requires its radial mode')
        return disk_profile(x,y,0.,speed=0.,amplitude=amplitude)
    s,f=mode
    x=x-L/2;y=y-L/2
    return .5+amplitude*np.real(np.exp(2j*np.arctan2(y,x))*sample(s,f,np.hypot(x,y)/R))


def qualify(folder):
    """Require separate spatial and temporal refinement before running OpenLB."""
    folder=Path(folder);folder.mkdir();cache={};comparisons=[];calls=[]
    probe=np.linspace(0.,1.,4097);limit=1e-4
    julia=runtime()
    def compute(n,steps):
        key=(n,steps)
        if key not in cache:
            p=folder/f'radial_{n}_{steps}.csv'
            cmd=[*julia,str(HERE/'radial_reference.jl'),str(p),*map(str,[n,steps,*TIMES,D,R,OMEGA0])]
            with p.with_suffix('.log').open('w') as stream:
                subprocess.run(cmd,check=True,stdout=stream,stderr=subprocess.STDOUT,timeout=120)
            calls.append(cmd);s,f=read_mode(p)
            if len(s)!=n or not np.isfinite(f).all() or np.max(abs(f))>1+limit:raise ValueError('Invalid radial reference')
            cache[key]=p,s,f
        return cache[key]
    def compare(a,b,axis):
        _,x,f=compute(*a);_,y,g=compute(*b)
        errors=[float(np.max(abs(sample(x,fi,probe)-sample(y,gi,probe)))) for fi,gi in zip(f,g)]
        comparisons.append(dict(coarse=a,fine=b,axis=axis,maximum_over_amplitude=errors))
        return max(errors)<=limit
    for n in (512,1024,2048,4096):
        if not compare((n//2,4096),(n,4096),'radial_spacing'):continue
        for steps in (4096,8192):
            if not compare((n,steps//2),(n,steps),'time_step'):continue
            if steps!=4096 and not compare((n//2,steps),(n,steps),'radial_spacing_at_selected_time_step'):break
            path,s,f=compute(n,steps)
            metadata=dict(passed=True,reference_file=path.name,reference_sha256=sha(path),radial_points=n,
                steps_per_window=steps,tolerance=limit,times_s=TIMES,probe_radii=len(probe),comparisons=comparisons,
                julia=subprocess.check_output([*julia,'-e','print(VERSION)'],text=True,timeout=30),commands=calls)
            write_json(folder/'qualification.json',metadata)
            return s,f,metadata
    write_json(folder/'qualification.json',dict(passed=False,comparisons=comparisons,commands=calls))
    raise ValueError(f'Radial reference convergence failed: {comparisons}')
