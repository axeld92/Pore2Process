"""Actual stock-DWSIM pump with the low-head liquid route, in a private adapter."""
import json
from pathlib import Path
import subprocess
from pore2process_demo.runtime.dwsim import DOTNET, DISTRIBUTION, environment
from pore2process_demo.runtime.openlb import ROOT
from pore2process_demo.io import write_json
from pore2process_demo.io import sha

ARTIFACTS=ROOT/'.local/adsorption/dwsim_runner'
RUNNER=ARTIFACTS/'bin/AdsorptionPumpRunner/debug/AdsorptionPumpRunner.dll'

def build_adapter(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    sdk=subprocess.check_output([str(DOTNET),'--version'],env=environment(),text=True,timeout=30).strip()
    if sdk!='10.0.403':raise ValueError('Review changed .NET SDK: '+sdk)
    command=[str(DOTNET),'build',str(ROOT/'examples/adsorption/dwsim_runner/AdsorptionPumpRunner.csproj'),
        '--artifacts-path',str(ARTIFACTS),'--ignore-failed-sources','--disable-build-servers','-p:UseSharedCompilation=false']
    with (directory/'build.log').open('w') as f:subprocess.run(command,env=environment(),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
    return {'command':command,'dotnet_sdk':sdk,'binary_sha256':sha(RUNNER)}

def solve_case(directory,inputs):
    directory=Path(directory);directory.mkdir(parents=True)
    write_json(directory/'inputs.json',inputs)
    with (directory/'stdout.log').open('w') as f:
        subprocess.run([str(DOTNET),str(RUNNER),str(DISTRIBUTION),str((directory/'inputs.json').resolve()),str(directory.resolve())],
            env=environment(),cwd=DISTRIBUTION,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=180)
    value=json.loads((directory/'summary.json').read_text())
    if value['inputs']!=inputs:raise ValueError('Actual DWSIM input snapshot changed.')
    return value
