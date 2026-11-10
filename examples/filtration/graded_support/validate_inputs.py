"""Exercise native input and unconverged-run guards without long solves."""
import json
from pathlib import Path
import subprocess
from datetime import datetime,timezone

from pore2process_demo.io import sha
from examples.filtration.graded_support.build import BINARY
from examples.filtration.graded_support.model import ROOT
from pore2process_demo.io import write_json

def main():
    base=ROOT/"results/graded_support/input_validation"/datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    base.mkdir(parents=True)
    bad=(('bulk_too_small',['--bulk-radius-m','3e-6']),('bulk_too_large',['--bulk-radius-m','7e-6']),
         ('bulk_nan',['--bulk-radius-m','nan']),('bad_grid',['--resolution','17']),
         ('small_lead',['--lead-length-m','8e-6']),('nan_lead',['--lead-length-m','nan']),
         ('unaligned_lead',['--lead-length-m','16.1e-6']),('bad_neck',['--radius-m','5e-6']),
         ('zero_pressure',['--pressure-drop-Pa','0']),('large_pressure',['--pressure-drop-Pa','21']),
         ('invalid_tau',['--tau','.5']),('nonfinite_pressure',['--pressure-drop-Pa','inf']),
         ('invalid_viscosity',['--viscosity-Pa-s','0']),('invalid_density',['--density-kg-m3','-1']),
         ('invalid_tolerance',['--steady-tolerance','nan']),('trailing_numeric',['--bulk-radius-m','5e-6x']),
         ('missing_value',['--bulk-radius-m']),('unconverged',['--resolution','8','--max-steps','1']))
    cases=[]
    for name,options in bad:
        directory=base/name;directory.mkdir()
        command=[str(BINARY),*options]
        result=subprocess.run(command,cwd=directory,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=20)
        (directory/'stdout.log').write_text(result.stdout)
        if result.returncode==0 or (directory/'raw.json').exists() or 'Graded support failed:' not in result.stdout:
            raise RuntimeError(f'Native guard did not reject cleanly: {name}')
        if name=='unconverged' and 'Maximum steps reached' not in result.stdout:
            raise RuntimeError('Expected an actual nonconvergence rejection.')
        cases.append({'name':name,'command':command,'returncode':result.returncode,'no_successful_raw_export':True})
    summary={'status':'passed','cases':cases,'binary_sha256':sha(BINARY),'validator_sha256':sha(Path(__file__)),
        'scope':'Actual native executable rejected 17 invalid inputs and one deliberately unconverged case; no successful simulation claimed.'}
    write_json(base/'summary.json',summary)
    print(f'{len(cases)} native rejection checks passed: {base}')

if __name__=='__main__': main()
