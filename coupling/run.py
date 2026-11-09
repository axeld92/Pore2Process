"""Real OpenLB → Python/Julia → DWSIM hydraulic design comparison.

Run from the project root: uv run python -m coupling.run
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter
from uuid import uuid4

import numpy as np

from compatibility.run_dwsim import build_adapter, solve_case, DISTRIBUTION, ARTIFACTS
from hydraulic.build import verify_vendor
from hydraulic.run import simulate, analyze, write_json, ROOT, BINARY

# Constant-property, low-Re water demonstration; no real material is claimed.
PROCESS = {"temperature_K":298.15,"inlet_absolute_pressure_Pa":1e5,
           "column_length_m":1.,"column_area_m2":1.,"volume_flow_m3_s":1e-5,
           "pump_efficiency_fraction":0.75}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


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


def execute(output):
    run_id=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    directory=output / "runs" / run_id
    directory.mkdir(parents=True)
    build=build_adapter(directory / "adapter_build")
    # Ask the actual DWSIM water property package for the fluid anchor first.
    # This is a property probe, and is not one of the material-design results.
    anchor_inputs={"Temperature_K":PROCESS["temperature_K"],
        "Inlet_pressure_Pa":PROCESS["inlet_absolute_pressure_Pa"],"Outlet_pressure_Pa":101000.,
        "Mass_flow_kg_s":0.01,"Pump_efficiency_fraction":PROCESS["pump_efficiency_fraction"],
        "Scope":"Water property anchor for the connected hydraulic demonstration."}
    print("DWSIM: obtaining liquid-water properties at the shared inlet state",flush=True)
    anchor=solve_case(directory / "fluid_anchor",anchor_inputs)
    rho,mu=anchor["inlet_density_kg_m3"],anchor["inlet_viscosity_Pa_s"]
    if anchor["inlet_vapor_fraction"]!=0 or not np.isfinite(rho+mu) or min(rho,mu)<=0:
        raise RuntimeError("Property anchor is not finite, positive, liquid water.")
    from pore2process_demo.bridge import runtime
    started=perf_counter()
    jl=runtime()
    initialization_s=perf_counter()-started
    cases=[]
    for name,radius in (("radius_50um",50e-6),("radius_40um",40e-6)):
        pore=[]
        # Independent solves of both physical geometries at two resolutions.
        for n in (24,40):
            path=directory / name / f"n{n}"
            raw,command=simulate(path,n,radius=radius,viscosity=mu,density=rho)
            metrics=analyze(path,raw,jl)
            record={"resolution":n,"inputs_and_measurements":raw,"metrics":metrics,"command":command,
                    "artifacts_sha256":{p.name:sha(p) for p in path.iterdir() if p.suffix in (".csv",".json")}}
            write_json(path / "analysis.json",record)
            pore.append(record)
            print(f"  {name}, N={n}: k={metrics['permeability_cell_m2']:.9g} m²; error={metrics['permeability_relative_error']:.3%}",flush=True)
        fine=pore[-1]
        k=fine["metrics"]["permeability_cell_m2"]
        dp=float(jl.predict_column_pressure_drop(k,mu,PROCESS["column_length_m"],
                                                PROCESS["column_area_m2"],PROCESS["volume_flow_m3_s"]))
        porosity=fine["metrics"]["porosity"]
        # Independent estimate verifies that scale translation stays in low-Re flow.
        max_velocity=2*PROCESS["volume_flow_m3_s"]/(PROCESS["column_area_m2"]*porosity)
        max_re=rho*max_velocity*2*radius/mu
        inputs={"Temperature_K":PROCESS["temperature_K"],
                "Inlet_pressure_Pa":PROCESS["inlet_absolute_pressure_Pa"],
                "Outlet_pressure_Pa":PROCESS["inlet_absolute_pressure_Pa"]+dp,
                "Mass_flow_kg_s":rho*PROCESS["volume_flow_m3_s"],
                "Pump_efficiency_fraction":PROCESS["pump_efficiency_fraction"],
                "Scope":f"OpenLB-derived permeability, Julia Darcy load, fixed-flow pump demonstration; {run_id}/{name}."}
        translation={"permeability_m2":k,"viscosity_Pa_s":mu,"density_kg_m3":rho,
                     "pressure_drop_Pa":dp,"process":PROCESS,"pore_run":f"{run_id}/{name}/n40",
                     "estimated_max_pore_velocity_m_s":max_velocity,"estimated_max_pore_Re":max_re,
                     "assumptions":"Aligned repeating straight pores; fully saturated Newtonian water; constant inlet properties; no gravity, other losses or entrance loss. Pump pressure rise equals the separately calculated porous load; no porous unit operation in DWSIM."}
        write_json(directory / name / "translation.json",translation)
        print(f"DWSIM: solving pump for {name}; porous pressure requirement {dp:.6g} Pa",flush=True)
        process=solve_case(directory / name / "dwsim",inputs)
        checks=validate_process(process,inputs,rho,mu,dp)
        checks.update({"fine_permeability_error_below_1_percent":fine["metrics"]["permeability_relative_error"]<0.01,
                       "fine_profile_error_below_1_percent":fine["metrics"]["velocity_profile_relative_L2"]<0.01,
                       "refinement_reduces_permeability_error":fine["metrics"]["permeability_relative_error"]<pore[0]["metrics"]["permeability_relative_error"],
                       "fine_mass_spread_below_0_1_percent":fine["metrics"]["all_plane_mass_relative_spread"]<0.001,
                       "process_max_pore_Re_below_0_05":max_re<0.05})
        if not all(checks.values()):
            raise RuntimeError(f"Design acceptance failed: {name}, {checks}")
        case={"name":name,"radius_m":radius,"pore_cases":pore,"translation":translation,
              "dwsim":process,"checks":checks}
        write_json(directory / name / "summary.json",case)
        cases.append(case)
    first,second=cases
    kratio=second["translation"]["permeability_m2"]/first["translation"]["permeability_m2"]
    pratio=second["translation"]["pressure_drop_Pa"]/first["translation"]["pressure_drop_Pa"]
    wratio=second["dwsim"]["shaft_power_W"]/first["dwsim"]["shaft_power_W"]
    reference=(second["radius_m"]/first["radius_m"])**4
    checks={"fixed_representative_area":all(c["pore_cases"][-1]["metrics"]["representative_area_m2"]==4e-8 for c in cases),
            "permeability_ratio_matches_radius_fourth_power_within_1_percent":abs(kratio/reference-1)<0.01,
            "actual_power_ratio_equals_inverse_measured_k_ratio":np.isclose(wratio,1/kratio,rtol=1e-9),
            "pressure_ratio_equals_inverse_measured_k_ratio":np.isclose(pratio,1/kratio,rtol=1e-12),
            "smaller_pore_increases_actual_DWSIM_power":wratio>1,
            "all_design_checks":all(all(c["checks"].values()) for c in cases)}
    if not all(checks.values()):
        raise RuntimeError(f"End-to-end ratio checks failed: {checks}")
    names=("hydraulic/openlb/capillary3d.cpp","hydraulic/openlb/case.h","hydraulic/openlb/config.mk",
           "hydraulic/build.py","hydraulic/run.py","coupling/run.py","julia/src/Hydraulics.jl",
           "compatibility/run_dwsim.py","compatibility/dwsim_runner/Program.cs",
           "src/pore2process_demo/bridge.py","uv.lock","julia/Project.toml","julia/Manifest.toml")
    summary={"schema_version":1,"run_id":run_id,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
             "status":"passed","process_inputs":PROCESS,"fluid_anchor":anchor,"designs":cases,
             "ratios":{"permeability_40_to_50":kratio,"pressure_40_to_50":pratio,
                       "actual_DWSIM_power_40_to_50":wratio,"analytical_radius_fourth_power":reference},
             "checks":{k:bool(v) for k,v in checks.items()},"julia_initialization_s":initialization_s,
             "versions":{"python":sys.version.split()[0],"julia":str(jl.seval("VERSION")),
                         "openlb":"1.9r0","dwsim":first["dwsim"]["version"],**build},
             "openlb_build":json.loads((ROOT / "results/hydraulic/build.json").read_text()),
             "vendor_integrity":verify_vendor(),"binary_sha256":sha(BINARY),
             "implementation_sha256":{name:sha(ROOT / name) for name in names},
             "runtime_sha256":{str(p.relative_to(ROOT)):sha(p) for p in
                (ARTIFACTS / "bin/DwsimProbe/debug/DwsimProbe.dll",DISTRIBUTION / "DWSIM.UnitOperations.dll",DISTRIBUTION / "DWSIM.Thermodynamics.dll")},
             "contract":"All boundaries SI. OpenLB gauge pressure and +x flux; representative area fixed at 4e-8 m²; x=i*dx, y=(j+0.5)*dx,z=(k+0.5)*dx, zero-based indices. DWSIM absolute pressure=inlet absolute + Julia porous pressure loss; mass flow=rho*Q. Internal kW/kJ outputs converted to W/J."}
    # Hash every raw export, input, translation and saved flowsheet for replay checks.
    summary["artifacts_sha256"]={str(p.relative_to(directory)):sha(p) for p in directory.rglob("*")
        if p.is_file() and p.suffix in (".csv",".json",".dwxmz")}
    write_json(directory / "summary.json",summary)
    from coupling.plot import render
    render(directory / "summary.json")
    write_json(output / "latest.json",{"run_id":run_id,"summary":str(directory / "summary.json"),"status":"passed"})
    print(f"Connected workflow passed; power ratio={wratio:.7g}: {directory / 'summary.json'}",flush=True)
    return directory / "summary.json"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT / "results/coupling")
    args=parser.parse_args()
    # Build is intentionally external to the original OpenLB distribution.
    subprocess.run([sys.executable,str(ROOT / "hydraulic/build.py")],check=True,timeout=700)
    execute(args.output.resolve())


if __name__=="__main__":
    main()
