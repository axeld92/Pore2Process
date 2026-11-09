"""Build a private .NET adapter and solve/check actual DWSIM flowsheets."""
from datetime import datetime, timezone
from pathlib import Path
import json
import math
import os
import subprocess

ROOT = Path(__file__).resolve().parents[1]
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
    build = [str(DOTNET), "build", str(ROOT / "compatibility/dwsim_runner/DwsimProbe.csproj"),
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


def main():
    build_record=build_adapter(OUTPUT)
    cases = []
    for name, pressure in (("baseline", 1e6), ("double_pressure_rise", 1.9e6)):
        directory = OUTPUT / name
        directory.mkdir(parents=True, exist_ok=True)
        inputs = {"Temperature_K": 298.15, "Inlet_pressure_Pa": 1e5,
                  "Outlet_pressure_Pa": pressure, "Mass_flow_kg_s": 10.0,
                  "Pump_efficiency_fraction": 0.75,
                  "Scope": "Independent water-pump compatibility probe; no OpenLB-derived input yet."}
        print(f"Solving real DWSIM water-pump case: {name}...", flush=True)
        result=solve_case(directory,inputs)
        numeric = [value for value in result.values() if isinstance(value, (int, float))]
        if not all(math.isfinite(value) for value in numeric):
            raise RuntimeError("Non-finite DWSIM result.")
        checks = [abs(result["outlet_pressure_Pa"] - pressure) <= 1e-6,
                  abs(result["mass_balance_error_kg_s"]) <= 1e-9 * inputs["Mass_flow_kg_s"],
                  abs(result["energy_balance_error_W"]) <= 1e-6 * result["shaft_power_W"],
                  result["pump_power_relative_error"] <= 1e-9,
                  abs(result["energy_stream_power_W"] - result["shaft_power_W"]) <= 1e-6,
                  result["inlet_vapor_fraction"] == 0, result["outlet_vapor_fraction"] == 0,
                  result["save_reload_solve_checked"] is True,
                  (directory / "water_pump.dwxmz").is_file()]
        if not all(checks):
            raise RuntimeError(f"DWSIM acceptance failed for {name}: {checks}; inspect summary/log.")
        cases.append(result)
        print(f"{name}: {result['shaft_power_W']:.6f} W; checks passed", flush=True)
    ratio = cases[1]["shaft_power_W"] / cases[0]["shaft_power_W"]
    if not math.isclose(ratio, 2.0, rel_tol=1e-9):
        raise RuntimeError("Doubling pressure rise did not double pump power.")
    summary = {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "cases": cases,
               "power_ratio": ratio, **build_record,
               "contract": "JSON field suffixes give SI units; efficiency is a fraction; internal kW/kJ converted at adapter boundary"}
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(f"Passed: {OUTPUT / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
