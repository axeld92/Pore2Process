"""Measure a synthetic FCC bed; keep its geometry distinct from experimental carbon."""
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from hydraulic.build import ROOT,build_case,verify_vendor
from hydraulic.run import write_json
from coupling.run import sha

BUILD=ROOT/'.local/adsorption/bed'
def read_json(path):return json.loads(Path(path).read_text())

def analyze(directory):
    p=Path(directory);r=json.loads((p/'raw.json').read_text());n=r['resolution'];h=r['dx_m']
    if not r['converged']:raise ValueError('Unconverged native bed.')
    a=np.atleast_1d(np.genfromtxt(p/'field.csv',delimiter=',',names=True))
    if len(a)!=r['fluid_nodes'] or len(set(zip(a['i'],a['j'],a['k'])))!=len(a):raise ValueError('Native indexing/count mismatch.')
    u=float(np.sum(a['ux_m_s'])/n**3);K=r['viscosity_Pa_s']*u/r['gradient_Pa_m']
    flux=np.array([np.sum(a['density_kg_m3'][a['i']==i]*a['ux_m_s'][a['i']==i])*h*h for i in range(n)])
    hist=np.genfromtxt(p/'history.csv',delimiter=',',names=True)
    residuals={key:float(np.std(hist[key][-6:],ddof=1)/abs(np.mean(hist[key][-6:]))) for key in ('superficial_velocity_m_s','mean_lattice_energy')}
    return {'raw':r,'permeability_m2':K,'superficial_velocity_m_s':u,'voxel_porosity':len(a)/n**3,
        'particle_Re_superficial':r['density_kg_m3']*abs(u)*r['particle_diameter_m']/r['viscosity_Pa_s'],
        'plane_mass_relative_spread':float(np.ptp(flux)/abs(np.mean(flux))),'stationarity':residuals}

def main():
    base=ROOT/'results/adsorption/bed';base.mkdir(parents=True,exist_ok=True)
    build=build_case(ROOT/'adsorption/openlb_bed',BUILD,base/'build');cases={}
    tau=lambda n:.5+.0009*n
    jobs=[('grid24',24,10,tau(24),.0023),('grid32',32,10,tau(32),.0023),('grid48',48,10,tau(48),.0023),
          ('pressure',48,20,tau(48),.0023),('time',48,10,.525,.0023),('scaled',48,40,tau(48),.00115),
          ('smallest',48,160,tau(48),.000575),('cold',24,10,tau(24),.0023)]
    for name,n,g,tau,d in jobs:
        p=base/name;p.mkdir(exist_ok=True);command=[str(BUILD/'capillary3d'),str(n),str(g),str(tau),str(d)]
        print('Native FCC bed: '+name,flush=True)
        initial=ROOT/'results/adsorption/bed_initial_compressibility_rejected'/('grid'+str(n))/'field.csv'
        scale=g/10*(d/.0023)**2
        warm=initial.exists() and name!='cold'
        if warm:command += [str(initial),str(scale)]
        request={'command':command,'binary_sha256':sha(BUILD/'capillary3d'),'initial_field_sha256':sha(initial) if warm else None,
            'initial_field_scope':'Optional native velocity estimate from the rejected acoustic scaling; every new case must fully converge with revised dt and meet all unchanged screens.'}
        if (p/'request.json').exists() and read_json(p/'request.json')!=request:raise ValueError('Native request changed; use a fresh output folder.')
        if not (p/'raw.json').exists() or not read_json(p/'raw.json')['converged']:
            with (p/'stdout.log').open('w') as f:subprocess.run(command,cwd=p,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=1800)
        write_json(p/'request.json',request);cases[name]=analyze(p)
    K=lambda name:cases[name]['permeability_m2']
    checks={'grid':abs(K('grid48')/K('grid32')-1)<.05,'pressure':abs(K('pressure')/K('grid48')-1)<.02,
        'time':abs(K('time')/K('grid48')-1)<.02,'geometric_scale':abs(4*K('scaled')/K('grid48')-1)<.02,
        'smallest_scale':abs(16*K('smallest')/K('grid48')-1)<.02,
        'initial_state':abs(K('cold')/K('grid24')-1)<1e-5,
        'mass':all(c['plane_mass_relative_spread']<1e-6 for c in cases.values()),
        'Mach':all(c['raw']['max_Mach']<.05 for c in cases.values()),
        'density':all(c['raw']['max_density_departure']<1e-4 for c in cases.values()),
        'porosity':all(abs(c['voxel_porosity']-.6)<.015 for c in cases.values()),
        'stationarity':all(max(c['stationarity'].values())<=1.0001e-7 for c in cases.values())}
    if not all(checks.values()):raise ValueError('FCC bed verification failed: '+str(checks))
    # Ergun is a comparison, not the oracle for an ordered FCC arrangement.
    ergun=.0023**2*.6**3/(150*(1-.6)**2)
    s={'status':'passed','cases':cases,'accepted_case':'grid48','permeability_m2':K('grid48'),
        'particle_diameter_m':.0023,'declared_porosity':.6,'checks':checks,'build':build,'vendor':verify_vendor(),
        'Ergun_comparison_only_m2':ergun,'scope':'Periodic ordered FCC impermeable spheres; unresolved intraparticle adsorption is a separate Julia model. Not reconstruction or calibration of experimental carbon.'}
    write_json(base/'summary.json',s);print(json.dumps({'checks':checks,'K_m2':K('grid48')},indent=2))
if __name__=='__main__':main()
