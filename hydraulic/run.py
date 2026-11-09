"""Run the real 3D OpenLB suite; use Julia to infer permeability from measurements."""
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

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / ".local/hydraulic/build/capillary3d"


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def simulate(directory, resolution, pressure=0.1, tau=0.8, tolerance=1e-8,
             radius=50e-6, viscosity=1e-3, density=1000., *, binary=BINARY, options=(), timeout_s=600):
    directory.mkdir(parents=True)
    command = [str(binary), "--resolution", str(resolution), "--pressure-drop-Pa", str(pressure),
               "--tau", str(tau), "--steady-tolerance", str(tolerance), "--max-steps", "60000",
               "--radius-m", str(radius), "--viscosity-Pa-s", str(viscosity), "--density-kg-m3", str(density)]
    command.extend(options)
    print(f"OpenLB: {directory.name}, N={resolution}, Δp={pressure} Pa, tau={tau}", flush=True)
    with (directory / "stdout.log").open("w") as log:
        result = subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_s)
    if result.returncode:
        raise RuntimeError(f"OpenLB failed ({result.returncode}); inspect {directory / 'stdout.log'}")
    raw = json.loads((directory / "raw.json").read_text())
    if not raw["converged"] or raw["resolution"] != resolution or raw["imposed_delta_p_Pa"] != pressure:
        raise RuntimeError("Solver result does not match the requested input snapshot.")
    if raw["tau"] != tau or raw["steady_tolerance"] != tolerance:
        raise RuntimeError("Solver relaxation/convergence settings changed.")
    if (raw["diameter_m"],raw["viscosity_Pa_s"],raw["reference_density_kg_m3"]) != (2*radius,viscosity,density):
        raise RuntimeError("Solver geometry/fluid snapshot mismatch.")
    numeric = [value for value in raw.values() if isinstance(value,(int,float))]
    if not all(np.isfinite(numeric)):
        raise RuntimeError("Non-finite OpenLB output.")
    return raw, command


def analyze(directory, raw, jl):
    planes = np.genfromtxt(directory / "planes.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    samples = np.genfromtxt(directory / "midplane.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    mu, diameter, area = raw["viscosity_Pa_s"], raw["diameter_m"], raw["cell_width_m"]**2
    measured_k = float(jl.infer_permeability(raw["tap_length_m"], area, mu,
                                            raw["tap_delta_p_Pa"], raw["Q_mid_m3_s"]))
    reference = jl.capillary_reference(diameter/2, raw["length_m"], mu, raw["imposed_delta_p_Pa"], area)
    q_ref, k_ref = float(reference.flow_m3_s), float(reference.permeability_cell_m2)
    radius2 = (samples["y_m"]-diameter/2)**2 + (samples["z_m"]-diameter/2)**2
    gradient = raw["tap_delta_p_Pa"] / raw["tap_length_m"]
    u_ref = gradient/(4*mu) * ((diameter/2)**2-radius2)
    inner = planes[(planes["x_m"] >= raw["length_m"]/4-raw["dx_m"]*0.01) &
                   (planes["x_m"] <= raw["length_m"]*3/4+raw["dx_m"]*0.01)]
    mass = planes["mass_flow_kg_s"]
    mass_scale = np.mean(np.abs(mass))
    inner_mass = inner["mass_flow_kg_s"]
    p_line = inner["mean_pressure_Pa"][0] - gradient*(inner["x_m"]-inner["x_m"][0])
    metrics = {
        "permeability_cell_m2": measured_k, "representative_area_m2": area,
        "reference_Q_m3_s": q_ref, "reference_permeability_cell_m2": k_ref,
        "reference_permeability_pore_area_m2": float(reference.permeability_pore_m2),
        "porosity": float(reference.porosity),
        "Q_relative_error": abs(raw["Q_mid_m3_s"]/q_ref-1),
        "permeability_relative_error": abs(measured_k/k_ref-1),
        "velocity_profile_relative_L2": float(np.linalg.norm(samples["ux_m_s"]-u_ref)/np.linalg.norm(u_ref)),
        "transverse_velocity_ratio": float(np.max(np.hypot(samples["uy_m_s"],samples["uz_m_s"]))/np.max(np.abs(u_ref))),
        "pressure_linearity_relative_error": float(np.max(np.abs(inner["mean_pressure_Pa"]-p_line))/abs(raw["tap_delta_p_Pa"])),
        "inlet_outlet_mass_relative_error": abs(raw["mass_in_kg_s"]-raw["mass_out_kg_s"])/mass_scale,
        "all_plane_mass_relative_spread": float(np.ptp(mass)/mass_scale),
        "inner_plane_mass_relative_spread": float(np.ptp(inner_mass)/np.mean(np.abs(inner_mass))),
        "inlet_outlet_volume_relative_difference": abs(raw["Q_in_m3_s"]-raw["Q_out_m3_s"])/abs(raw["Q_mid_m3_s"]),
        "max_density_departure_fraction": float(np.max(np.abs(samples["density_kg_m3"]/raw["reference_density_kg_m3"]-1))),
        "measured_max_mach": float(np.max(np.abs(samples["ux_m_s"]))*raw["dt_s"]/raw["dx_m"]*np.sqrt(3)),
        "measured_max_velocity_Re": float(np.max(np.abs(samples["ux_m_s"]))*diameter*raw["reference_density_kg_m3"]/mu),
        # Keep the same physical pressure gradient as the pore benchmark.
        "darcy_column_flow_m3_s": float(jl.predict_column_flow(measured_k,mu,0.02,1e-4,5.0)),
    }
    if raw.get("convergence_method","").startswith("OpenLB ValueTracer"):
        history=np.atleast_1d(np.genfromtxt(directory / "history.csv",delimiter=",",names=True))
        window=raw["convergence_window_samples"]
        if len(history)<window:
            raise RuntimeError("Convergence history does not cover the declared OpenLB tracer window.")
        for field,label in (("Q_mid_m3_s","Q"),("delta_p_taps_Pa","dp")):
            values=history[field][-window:]
            relative_std=float(np.std(values,ddof=1)/abs(np.mean(values)))
            metrics[f"convergence_{label}_relative_std"]=relative_std
            if not np.isfinite(relative_std) or relative_std>=raw["steady_tolerance"]*(1+1e-5):
                raise RuntimeError("Exported history fails the OpenLB ValueTracer stop criterion.")
    if not all(np.isfinite(list(metrics.values()))):
        raise RuntimeError("Non-finite post-processing result.")
    nx=raw.get("axial_intervals",4*raw["resolution"])
    offset=raw.get("transverse_index_offset",0.)
    if len(planes) != nx+1 or len(samples) != raw["inlet_nodes"]:
        raise RuntimeError("Unexpected array shapes / boundary-to-midplane indexing.")
    if any(samples[name].dtype != np.dtype("int64") for name in ("i","j","k")):
        raise RuntimeError("CSV lattice coordinates must materialize as int64.")
    for index, coordinate in (("i","x_m"),("j","y_m"),("k","z_m")):
        shift=0. if index=="i" else offset
        if not np.allclose((samples[index]+shift)*raw["dx_m"],samples[coordinate],rtol=0,atol=1e-12*raw["dx_m"]):
            raise RuntimeError("Exported zero-based indices and physical coordinates disagree.")
    if not np.array_equal(planes["i"],np.arange(len(planes))) or not np.all(samples["i"] == nx//2):
        raise RuntimeError("Plane order or midplane indices violate the contract.")
    if abs(raw["inlet_pressure_Pa"]-raw["imposed_delta_p_Pa"]) > 1e-9 or abs(raw["outlet_pressure_Pa"]) > 1e-9:
        raise RuntimeError("Pressure boundary values differ from the input.")
    if metrics["inlet_outlet_mass_relative_error"] > 2e-5 or metrics["inner_plane_mass_relative_spread"] > 1e-4:
        raise RuntimeError(f"Mass-flow acceptance failed: {directory}; {metrics}")
    if metrics["all_plane_mass_relative_spread"] > 0.005 or metrics["pressure_linearity_relative_error"] > 0.002:
        raise RuntimeError(f"Plane-flow/pressure-profile acceptance failed: {directory}; {metrics}")
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/hydraulic")
    parser.add_argument("--reuse", type=Path, help="Reuse unchanged solver artifacts from this summary; binary, inputs and file hashes must match.")
    args = parser.parse_args()
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    output = args.output.resolve() / "runs" / run_id
    output.mkdir(parents=True)
    if not BINARY.is_file():
        raise SystemExit("Build the benchmark first: uv run python hydraulic/build.py")
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
               "implementation_sha256":{name:hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in (
                   "hydraulic/run.py","hydraulic/openlb/case.h","julia/src/Hydraulics.jl","src/pore2process_demo/bridge.py",
                   "uv.lock","julia/Project.toml","julia/Manifest.toml")},
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
