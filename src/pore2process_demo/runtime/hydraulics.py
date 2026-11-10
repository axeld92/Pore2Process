"""Existing native capillary measurements and hydraulic process checks."""
import json
import subprocess
from pathlib import Path
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import write_json
BINARY = ROOT / ".local/hydraulic/build/capillary3d"

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


def validate_process(result, inputs, rho, mu, dp):
    q=inputs["Mass_flow_kg_s"]/rho
    expected=q*dp/inputs["Pump_efficiency_fraction"]
    checks={
        "actual_input_snapshot":result["inputs"]==inputs,
        "absolute_outlet_pressure":abs(result["outlet_pressure_Pa"]-inputs["Outlet_pressure_Pa"])<1e-6,
        "liquid_water":result["inlet_vapor_fraction"]==0 and result["outlet_vapor_fraction"]==0,
        "same_density":np.isclose(result["inlet_density_kg_m3"],rho,rtol=1e-12),
        "same_viscosity":np.isclose(result["inlet_viscosity_Pa_s"],mu,rtol=1e-12),
        "constant_property_density_change_below_0_01_percent":abs(result["outlet_density_kg_m3"]/rho-1)<1e-4,
        "constant_property_viscosity_change_below_0_01_percent":abs(result["outlet_viscosity_Pa_s"]/mu-1)<1e-4,
        "mass_conservation":abs(result["mass_balance_error_kg_s"])<1e-10*inputs["Mass_flow_kg_s"],
        # Small pressure rises produce very small enthalpy increments. The flash
        # residual is reported explicitly; never equate it with exact energy closure.
        "energy_balance_below_0_1_percent_of_power":abs(result["energy_balance_error_W"])<1e-3*expected,
        "shaft_power_matches_hydraulics":np.isclose(result["shaft_power_W"],expected,rtol=1e-9),
        "energy_stream_matches_shaft_power":np.isclose(result["energy_stream_power_W"],expected,rtol=1e-9),
        "saved_reloaded_solved":result["save_reload_solve_checked"] is True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"DWSIM coupling acceptance failed: {checks}")
    return {name:bool(value) for name,value in checks.items()}
