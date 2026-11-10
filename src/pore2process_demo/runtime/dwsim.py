"""Build a private .NET adapter and solve/check actual DWSIM flowsheets."""
from datetime import datetime, timezone
from pathlib import Path
import json
import math
import os
import subprocess

from pore2process_demo.paths import ROOT
CACHE = ROOT / ".local/compatibility"
DOTNET = CACHE / "dotnet/dotnet"
DISTRIBUTION = CACHE / "dwsim"
OUTPUT = ROOT / "results/compatibility/dwsim"
ARTIFACTS = CACHE / "runner"


def environment():
    return {**os.environ, "DOTNET_ROOT": str(DOTNET.parent),
            "DOTNET_CLI_HOME": str(CACHE / "dotnet-home"),
            "DOTNET_CLI_TELEMETRY_OPTOUT": "1", "DOTNET_NOLOGO": "1",
            "DOTNET_GENERATE_ASPNET_CERTIFICATE": "false",
            "NUGET_PACKAGES": str(CACHE / "nuget"),
            "DYLD_LIBRARY_PATH": str(DISTRIBUTION)}


def build_adapter(output):
    output.mkdir(parents=True, exist_ok=True)
    env = environment()
    sdk_version = subprocess.check_output([str(DOTNET), "--version"], env=env, text=True, timeout=15).strip()
    if sdk_version != "10.0.403":
        raise SystemExit(f"Expected private SDK 10.0.403; found {sdk_version}. Recheck bootstrap.")
    build = [str(DOTNET), "build", str(ROOT / "src/pore2process_demo/runtime/dwsim_runner/DwsimProbe.csproj"),
             "--artifacts-path", str(ARTIFACTS), "--ignore-failed-sources",
             "--disable-build-servers", "-p:UseSharedCompilation=false"]
    with (output / "build.log").open("w") as log:
        subprocess.run(build, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
    return {"build_command":build,"dotnet_sdk":sdk_version}


def solve_case(directory, inputs):
    directory.mkdir(parents=True, exist_ok=True)
    input_file=directory / "inputs.json"
    input_file.write_text(json.dumps(inputs,indent=2,allow_nan=False)+"\n")
    command=[str(DOTNET),str(ARTIFACTS / "bin/DwsimProbe/debug/DwsimProbe.dll"),
             str(DISTRIBUTION),str(input_file.resolve()),str(directory.resolve())]
    with (directory / "stdout.log").open("w") as log:
        subprocess.run(command,env=environment(),cwd=DISTRIBUTION,stdout=log,
                       stderr=subprocess.STDOUT,check=True,timeout=180)
    result=json.loads((directory / "summary.json").read_text())
    if result["inputs"] != inputs:
        raise RuntimeError("DWSIM input snapshot mismatch.")
    return result
