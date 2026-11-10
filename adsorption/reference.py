"""Publish a compact real-simulator fixture, leaving full time series in results."""
import shutil
from pathlib import Path
from hydraulic.build import ROOT
from hydraulic.run import write_json
from .run import read
from coupling.run import sha

def main():
    p=Path(read(ROOT/'results/adsorption/process/latest.json')['summary']);s=read(p);v=read(p.parent/'verification.json')
    if s['status']!='passed_synthetic_example' or v['status']!='passed':raise ValueError('Verify the run before publishing its reference.')
    fixture={'summary_path':str(p),'summary_sha256':sha(p),'verification':v,'comparison':s['comparison'],
        'evidence':s['evidence'],'refinement':{k:s['refinement'][k] for k in ('space_breakthrough_relative','time_breakthrough_relative')},
        'cases':[{**{k:c[k] for k in ('particle_diameter_m','native_case','permeability_m2','inputs','dwsim')},
            'column_metrics':{k:v for k,v in c['column'].items() if not isinstance(v,list) and k!='total_s'}} for c in s['cases']]}
    write_json(ROOT/'docs/adsorption_reference.json',fixture);shutil.copy2(p.parent/'adsorption.png',ROOT/'docs/adsorption.png')
    print('Published actual adsorption-loop reference and figure.')
if __name__=='__main__':main()
