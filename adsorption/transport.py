"""Independent SI/Fourier check of stock OpenLB scalar TRT, before bed use."""
import json
from pathlib import Path
import subprocess
import numpy as np
from hydraulic.build import ROOT,build_case,verify_vendor
from coupling.run import sha
from hydraulic.run import write_json

BUILD=ROOT/'.local/adsorption/transport'

def analyze(directory):
    p=Path(directory);r=json.loads((p/'raw.json').read_text())
    a=np.genfromtxt(p/'field.csv',delimiter=',',names=True)
    t=r['steps']*r['dt_s'];k=2*np.pi/r['length_m'];x=a['x_m']
    exact=1+r['initial_amplitude']*np.sin(k*(x-r['velocity_m_s']*t))*np.exp(-r['diffusivity_m2_s']*k*k*t)
    scale=np.linalg.norm(exact-1) if r['initial_amplitude'] else np.linalg.norm(exact)
    relative=float(np.linalg.norm(a['normalized_concentration']-exact)/scale)
    amplitude=2*abs(np.mean((a['normalized_concentration']-1)*np.exp(-1j*k*x)))
    inferred=None if r['initial_amplitude']==0 else float(-np.log(amplitude/r['initial_amplitude'])/(k*k*t))
    return {'raw':r,'relative_profile_error':relative,'mean_concentration_error':abs(float(np.mean(a['normalized_concentration']))-1),
        'measured_diffusivity_m2_s':inferred,'minimum_concentration':float(min(a['normalized_concentration'])),
        'scope':'Prescribed-velocity periodic free-fluid verification; not a packed-bed dispersion or adsorption measurement.'}

def main():
    base=ROOT/'results/adsorption/transport';base.mkdir(parents=True,exist_ok=True)
    build=build_case(ROOT/'adsorption/openlb',BUILD,base/'build')
    cases={}
    for name,n,pe,tau,amp in [('grid16',16,4,.8,.2),('grid32',32,4,.8,.2),('grid64',64,4,.8,.2),
            ('reverse',64,-4,.8,.2),('diffusion',64,0,.8,.2),('time',64,4,.65,.2),('uniform',32,4,.8,0.)]:
        p=base/name;p.mkdir(exist_ok=True);command=[str(BUILD/'capillary3d'),str(n),str(pe),str(tau),str(amp)]
        with (p/'stdout.log').open('w') as f:subprocess.run(command,cwd=p,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
        write_json(p/'request.json',{'command':command,'binary_sha256':sha(BUILD/'capillary3d')})
        cases[name]=analyze(p)
    errors=[cases['grid'+str(n)]['relative_profile_error'] for n in (16,32,64)]
    orders=[float(np.log(errors[i]/errors[i+1])/np.log(2)) for i in range(2)]
    checks={'mass':all(c['mean_concentration_error']<1e-10 for c in cases.values()),
        'positive':all(c['minimum_concentration']>0 for c in cases.values()),
        'profile':all(c['relative_profile_error']<.005 for name,c in cases.items() if name!='grid16'),
        'refinement':min(orders)>1.7,'diffusivity':all(abs(c['measured_diffusivity_m2_s']/c['raw']['diffusivity_m2_s']-1)<.005 for name,c in cases.items() if c['measured_diffusivity_m2_s'] is not None and name!='grid16')}
    if not all(checks.values()):raise ValueError('Transport verification failed: '+str(checks)+' errors='+str(errors)+' orders='+str(orders))
    s={'status':'passed','cases':cases,'orders':orders,'checks':checks,'build':build,'vendor':verify_vendor()}
    write_json(base/'summary.json',s);print(json.dumps({'checks':checks,'orders':orders,'cases':len(cases)},indent=2))
if __name__=='__main__':main()
