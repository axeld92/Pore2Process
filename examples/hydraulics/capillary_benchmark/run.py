"""Run the real 3D OpenLB suite; use Julia to infer permeability from measurements."""
from pore2process_demo.provenance import implementation_files, source_manifest
from pore2process_demo.runtime.hydraulics import simulate, analyze
from pore2process_demo.io import write_json
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
BINARY = ROOT / ".local/hydraulic/build/capillary3d"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/hydraulic")
    parser.add_argument("--reuse", type=Path, help="Reuse unchanged solver artifacts from this summary; binary, inputs and file hashes must match.")
    args = parser.parse_args()
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    output = args.output.resolve() / "runs" / run_id
    output.mkdir(parents=True)
    if not BINARY.is_file():
        raise SystemExit("Build the benchmark first: uv run python -m examples.hydraulics.capillary_benchmark.build")
    binary_hash=hashlib.sha256(BINARY.read_bytes()).hexdigest()
    reuse_cases={}
    if args.reuse:
        previous=json.loads(args.reuse.read_text())
        if previous["binary_sha256"] != binary_hash:
            raise SystemExit("Cannot reuse results from a different solver binary.")
        reuse_cases={case["name"]:case for case in previous["cases"]}
    from pore2process_demo.bridge import runtime
    started = perf_counter()
    jl = runtime()
    initialization_s = perf_counter()-started
    cases = []
    for name, n, dp, tau, tolerance in (
        ("n16",16,0.1,0.8,1e-8), ("n24",24,0.1,0.8,1e-8), ("n40",40,0.1,0.8,1e-8),
        ("double_pressure",24,0.2,0.8,1e-8), ("reverse_pressure",24,-0.1,0.8,1e-8),
        ("half_time_step",24,0.1,0.65,1e-8), ("tighter_steady",24,0.1,0.8,1e-10),
        ("n32",32,0.1,0.8,1e-8),
    ):
        directory = output / name
        reused=None
        if name in reuse_cases:
            cached=reuse_cases[name]
            settings=cached["inputs_and_measurements"]
            expected=(n,dp,tau,tolerance)
            observed=tuple(settings[key] for key in ("resolution","imposed_delta_p_Pa","tau","steady_tolerance"))
            if observed != expected:
                raise RuntimeError(f"Cached input mismatch for {name}.")
            directory.mkdir()
            origin=args.reuse.resolve().parent / name
            for filename,digest in cached["artifacts_sha256"].items():
                source=origin / filename
                if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
                    raise RuntimeError(f"Cached solver artifact has changed: {source}")
                shutil.copy2(source,directory / filename)
            shutil.copy2(origin / "stdout.log",directory / "stdout.log")
            raw=json.loads((directory / "raw.json").read_text())
            command=cached["command"]
            reused={"run_id":cached["run_id"],"directory":str(origin)}
            print(f"Rechecking unchanged OpenLB artifacts: {name}",flush=True)
        else:
            raw, command = simulate(directory,n,dp,tau,tolerance,timeout_s=1800)
        metrics = analyze(directory,raw,jl)
        record = {"name":name,"run_id":run_id,"inputs_and_measurements":raw,"metrics":metrics,
                  "command":command,"reused_solver_artifacts":reused,
                  "artifacts_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in directory.iterdir() if p.suffix in (".csv",".json")}}
        write_json(directory / "summary.json",record)
        cases.append(record)
        print(f"  k={metrics['permeability_cell_m2']:.9g} m²; reference error={metrics['permeability_relative_error']:.3%}; mass checks passed",flush=True)
    base = cases[1]["metrics"]["permeability_cell_m2"]
    checks = {
        # Circular-wall discretization/quadrature errors can cancel on a coarse
        # grid. Report every grid; assess fine-grid stability and analytical error
        # instead of claiming monotonic error or a formal convergence order.
        "n32_to_n40_k_change_below_0_1_percent": abs(cases[2]["metrics"]["permeability_cell_m2"]/cases[7]["metrics"]["permeability_cell_m2"]-1) < 0.001,
        "n24_to_n40_k_change_below_0_5_percent": abs(cases[2]["metrics"]["permeability_cell_m2"]/base-1) < 0.005,
        "finest_k_error_below_1_percent": cases[2]["metrics"]["permeability_relative_error"] < 0.01,
        "finest_profile_error_below_1_percent": cases[2]["metrics"]["velocity_profile_relative_L2"] < 0.01,
        "finest_all_plane_mass_spread_below_0_1_percent": cases[2]["metrics"]["all_plane_mass_relative_spread"] < 0.001,
        "double_pressure_k_within_0_1_percent": abs(cases[3]["metrics"]["permeability_cell_m2"]/base-1) < 0.001,
        "reversal_k_within_0_1_percent": abs(cases[4]["metrics"]["permeability_cell_m2"]/base-1) < 0.001,
        "half_dt_k_within_0_5_percent": abs(cases[5]["metrics"]["permeability_cell_m2"]/base-1) < 0.005,
        "steady_tolerance_k_within_1e_7": abs(cases[6]["metrics"]["permeability_cell_m2"]/base-1) < 1e-7,
    }
    geometry_checks = []
    for n in (20,32,40):
        directory=output / f"geometry_n{n}"
        directory.mkdir()
        with (directory / "stdout.log").open("w") as log:
            subprocess.run([str(BINARY),"--resolution",str(n),"--geometry-only"],cwd=directory,
                           stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
        geometry_checks.append(json.loads((directory / "geometry.json").read_text()))
    summary = {"schema_version":1,"run_id":run_id,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
               "cases":cases,"checks":checks,"geometry_checks":geometry_checks,
               "julia_version":str(jl.seval("VERSION")),"python_version":sys.version.split()[0],
               "julia_initialization_s":initialization_s,
               "build":json.loads((ROOT / "results/hydraulic/build.json").read_text()),
               "binary_sha256":binary_hash,
               "implementation_sha256":source_manifest((
                   "examples/hydraulics/capillary_benchmark/run.py","examples/hydraulics/capillary_benchmark/openlb/case.h","julia/src/Hydraulics.jl","src/pore2process_demo/bridge.py",
                   "uv.lock","julia/Project.toml","julia/Manifest.toml")),
               "contract":"+x positive volume/mass flux on every plane; pressure is gauge Pa; representative area=cell_width_m²; CSV lattice indices zero-based; SI fields explicitly suffixed",
               "column_inputs":{"length_m":0.02,"area_m2":1e-4,"delta_p_Pa":5.0,"viscosity_Pa_s":0.001},
               "status":"passed" if all(checks.values()) else "failed"}
    write_json(output / "summary.json",summary)
    write_json(args.output.resolve() / "latest.json",{"run_id":run_id,"summary":str(output / "summary.json"),"status":summary["status"]})
    if not all(checks.values()):
        raise RuntimeError(f"Hydraulic suite failed: {checks}; inspect {output}")
    print(f"All hydraulic acceptance checks passed: {output / 'summary.json'}",flush=True)


if __name__ == "__main__":
    main()
