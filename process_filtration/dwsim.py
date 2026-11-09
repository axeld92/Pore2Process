"""Separate adapter; reuse the installed private SDK and unmodified DWSIM."""
import json
import hashlib
from pathlib import Path
import subprocess

from compatibility.run_dwsim import DOTNET, DISTRIBUTION, environment
from hydraulic.run import write_json

ROOT=Path(__file__).resolve().parents[1]
ARTIFACTS=ROOT/".local/process_filtration/runner"
RUNNER=ARTIFACTS/"bin/FiltrationRunner/debug/FiltrationRunner.dll"

def invoke(directory,input_path):
    command=[str(DOTNET),str(RUNNER),str(DISTRIBUTION),str(input_path),str(directory.resolve())]
    with (directory/"stdout.log").open("w") as f:
        subprocess.run(command,env=environment(),cwd=DISTRIBUTION,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=180)
    return command

def build(directory):
    directory.mkdir(parents=True,exist_ok=True)
    sdk=subprocess.check_output([str(DOTNET),"--version"],env=environment(),text=True,timeout=15).strip()
    if sdk!="10.0.403": raise RuntimeError(f"Review changed private SDK: {sdk}")
    command=[str(DOTNET),"build",str(ROOT/"process_filtration/dwsim_runner/FiltrationRunner.csproj"),
        "--artifacts-path",str(ARTIFACTS),"--ignore-failed-sources","--disable-build-servers","-p:UseSharedCompilation=false"]
    with (directory/"build.log").open("w") as f:
        subprocess.run(command,env=environment(),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
    invoke(directory,"--inspect")
    return {"build_command":command,"dotnet_sdk":sdk,"runner_sha256":hashlib.sha256(RUNNER.read_bytes()).hexdigest()}

def solve(directory,inputs):
    directory.mkdir(parents=True)
    input_path=(directory/"inputs.json").resolve();write_json(input_path,inputs)
    invoke(directory,input_path)
    result=json.loads((directory/"summary.json").read_text())
    if result["inputs"]!=inputs: raise RuntimeError("DWSIM input snapshot differs")
    return result
