"""Build and launch a separate stock-component DWSIM crossflow adapter."""
import hashlib
import json
from pathlib import Path
import subprocess

from pore2process_demo.runtime.dwsim import DOTNET, DISTRIBUTION, environment
from pore2process_demo.io import write_json

from pore2process_demo.paths import ROOT
ARTIFACTS=ROOT/'.local/crossflow_recycle/runner'
RUNNER=ARTIFACTS/'bin/CrossflowRunner/debug/CrossflowRunner.dll'

def invoke(directory,request,*,check=True):
    command=[str(DOTNET),str(RUNNER),str(DISTRIBUTION),str(request),str(directory.resolve())]
    with (directory/'stdout.log').open('w') as f:
        result=subprocess.run(command,env=environment(),cwd=DISTRIBUTION,stdout=f,stderr=subprocess.STDOUT,
                              timeout=180,check=check)
    return result.returncode

def build(directory):
    directory.mkdir(parents=True,exist_ok=True)
    sdk=subprocess.check_output([str(DOTNET),'--version'],env=environment(),text=True,timeout=15).strip()
    if sdk!='10.0.403':raise RuntimeError(f'Review changed SDK {sdk}')
    command=[str(DOTNET),'build',str(ROOT/'examples/filtration/crossflow_recycle/dwsim_runner/CrossflowRunner.csproj'),
        '--artifacts-path',str(ARTIFACTS),'--ignore-failed-sources','--disable-build-servers','-p:UseSharedCompilation=false']
    with (directory/'build.log').open('w') as f:
        subprocess.run(command,env=environment(),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
    invoke(directory,'--inspect')
    return {'build_command':command,'dotnet_sdk':sdk,'runner_sha256':hashlib.sha256(RUNNER.read_bytes()).hexdigest()}

def solve(directory,inputs):
    directory.mkdir(parents=True)
    request=(directory/'inputs.json').resolve();write_json(request,inputs)
    invoke(directory,request)
    result=json.loads((directory/'summary.json').read_text())
    if result['inputs']!=inputs:raise RuntimeError('DWSIM input snapshot differs')
    return result

def reject_nonconvergence(directory,inputs):
    directory.mkdir(parents=True)
    request=(directory/'inputs.json').resolve();write_json(request,inputs)
    code=invoke(directory,request,check=False)
    text=(directory/'stdout.log').read_text()
    if code==0 or 'RECYCLE_NOT_CONVERGED' not in text or (directory/'summary.json').exists():
        raise RuntimeError('Deliberately unconverged recycle was not rejected')
    record={'nonconvergence_rejected':True,'return_code':code,'maximum_iterations':inputs['Maximum_iterations']}
    write_json(directory/'rejection.json',record)
    return record
