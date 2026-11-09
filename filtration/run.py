"""Actual OpenLB supports → Julia fixed-flux layers → saved DWSIM pumps.

This is a synthetic clean-water sensitivity study, without experimental calibration.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from time import perf_counter
from uuid import uuid4

import numpy as np
from compatibility.run_dwsim import build_adapter, solve_case, ARTIFACTS, DISTRIBUTION
from coupling.run import sha, validate_process
from hydraulic.build import ROOT, verify_vendor
from hydraulic.run import simulate, analyze, write_json
from filtration.build import main as build_support, BINARY
from filtration.model import layer_load

SETTINGS = {"temperature_K":298.15, "inlet_absolute_pressure_Pa":1e5,
    "support_thickness_m":1e-3, "membrane_area_m2":1., "flux_LMH":100.,
    "flux_m_s":100.*1e-3/3600, "pump_efficiency_fraction":0.75,
    "other_fixed_resistance_m_inv":0.}
SCENARIOS = (("support_dominated",0.25),("balanced_baseline",1.),
             ("selective_dominated",10.),("strongly_selective_dominated",100.))


def pore_case(directory, n, radius, mu, rho, jl, pressure=10., length=4e-5):
    raw, command = simulate(directory,n,pressure=pressure,radius=radius,
        viscosity=mu,density=rho,binary=BINARY,options=("--length-m",str(length)))
    if not np.isclose(raw["length_m"],length,rtol=1e-14) or raw["cell_width_m"]!=2e-5:
        raise RuntimeError("Support geometry differs from its requested snapshot.")
    metrics=analyze(directory,raw,jl)
    checks={"k_error_below_1_percent":metrics["permeability_relative_error"]<0.01,
        "profile_error_below_1_percent":metrics["velocity_profile_relative_L2"]<0.01,
        "mass_spread_below_0_1_percent":metrics["all_plane_mass_relative_spread"]<0.001,
        "Re_below_0_05":metrics["measured_max_velocity_Re"]<0.05,
        "Mach_below_0_001":metrics["measured_max_mach"]<0.001}
    if not all(checks.values()):
        raise RuntimeError(f"Support acceptance failed: {directory}, {checks}")
    record={"inputs_and_measurements":raw,"metrics":metrics,"command":command,
            "checks":{k:bool(v) for k,v in checks.items()}}
    write_json(directory / "analysis.json",record)
    return record


def execute(output, reuse=None):
    run_id=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    directory=output / "runs" / run_id
    directory.mkdir(parents=True)
    adapter=build_adapter(directory / "adapter_build")
    anchor_inputs={"Temperature_K":SETTINGS["temperature_K"],"Inlet_pressure_Pa":1e5,
        "Outlet_pressure_Pa":101000.,"Mass_flow_kg_s":.01,"Pump_efficiency_fraction":.75,
        "Scope":"Liquid-water property anchor for the synthetic fixed-flux filtration study."}
    anchor=solve_case(directory / "fluid_anchor",anchor_inputs)
    rho,mu=anchor["inlet_density_kg_m3"],anchor["inlet_viscosity_Pa_s"]
    if anchor["inlet_vapor_fraction"]!=0 or not np.isfinite(rho+mu) or min(rho,mu)<=0:
        raise RuntimeError("The shared fluid anchor must be finite, positive liquid water.")
    from pore2process_demo.bridge import runtime
    start=perf_counter(); jl=runtime(); cold=perf_counter()-start
    # Optional reuse is only of a fully verified previous filtration run.
    previous=None
    if reuse:
        from filtration.check import verify
        verify(reuse)
        previous=json.loads(reuse.read_text())
        if previous["binary_sha256"]!=sha(BINARY) or previous["settings"]!=SETTINGS:
            raise RuntimeError("Cannot reuse a different binary or study settings.")
        if (previous["fluid_anchor"]["inlet_density_kg_m3"],previous["fluid_anchor"]["inlet_viscosity_Pa_s"])!=(rho,mu):
            raise RuntimeError("Fluid anchor changed; rerun pore simulations.")
    supports=[]
    for name,radius in (("baseline_4um",4e-6),("open_5um",5e-6)):
        cases=[]
        for n in (24,40):
            target=directory / "supports" / name / f"n{n}"
            if previous:
                import shutil
                shutil.copytree(reuse.parent / "supports" / name / f"n{n}",target)
                record=json.loads((target / "analysis.json").read_text())
            else:
                record=pore_case(target,n,radius,mu,rho,jl)
            cases.append(record)
            print(f"{name} N{n}: measured k={record['metrics']['permeability_cell_m2']:.10g} m²",flush=True)
        k=cases[-1]["metrics"]["permeability_cell_m2"]
        change=abs(k/cases[0]["metrics"]["permeability_cell_m2"]-1)
        if change>=.005:
            raise RuntimeError("N24→N40 support-k change exceeds 0.5%.")
        porosity=cases[-1]["metrics"]["porosity"]
        re=2*SETTINGS["flux_m_s"]/porosity*2*radius*rho/mu
        if re>=.05:
            raise RuntimeError("Upscaled support leaves the validated low-Re operating range.")
        supports.append({"name":name,"radius_m":radius,"permeability_m2":k,
            "porosity":porosity,"estimated_process_max_pore_Re":re,
            "grid_relative_change":change,"pore_cases":cases})
    extra=[]
    for name,pressure,length in (("low_pressure",1.,4e-5),("double_pressure",20.,4e-5),("double_domain",10.,8e-5)):
        target=directory / "validation" / name
        if previous and (reuse.parent / "validation" / name).exists():
            import shutil
            shutil.copytree(reuse.parent / "validation" / name,target)
            record=json.loads((target / "analysis.json").read_text())
        else:
            record=pore_case(target,24,4e-6,mu,rho,jl,pressure,length)
        base=supports[0]["pore_cases"][0]
        if abs(record["metrics"]["permeability_cell_m2"]/base["metrics"]["permeability_cell_m2"]-1)>=.001:
            raise RuntimeError(f"{name} changes inferred permeability by >=0.1%.")
        qratio=record["inputs_and_measurements"]["Q_mid_m3_s"]/base["inputs_and_measurements"]["Q_mid_m3_s"]
        expected=pressure/10.*4e-5/length
        if abs(qratio/expected-1)>=.001:
            raise RuntimeError(f"{name} flow scaling differs by >=0.1%.")
        extra.append({"name":name,**record,"Q_ratio":qratio})
    Rs0=SETTINGS["support_thickness_m"]/supports[0]["permeability_m2"]
    rows=[]; scenarios=[]
    for name,alpha in SCENARIOS:
        Rl=alpha*Rs0
        cases=[]
        for support in supports:
            load=layer_load(support["permeability_m2"],SETTINGS,Rl,mu)
            inputs={"Temperature_K":SETTINGS["temperature_K"],"Inlet_pressure_Pa":1e5,
                "Outlet_pressure_Pa":1e5+load["total_pressure_drop_Pa"],
                "Mass_flow_kg_s":rho*load["flow_m3_s"],"Pump_efficiency_fraction":.75,
                "Scope":f"Synthetic fixed-flux filtration: {run_id}/{name}/{support['name']}; externally calculated support + unchanged selective-layer load. No membrane separation unit."}
            target=directory / "scenarios" / name / support["name"]
            print(f"DWSIM {name}/{support['name']}: Δp={load['total_pressure_drop_Pa']:.6g} Pa",flush=True)
            process=solve_case(target,inputs)
            checks=validate_process(process,inputs,rho,mu,load["total_pressure_drop_Pa"])
            if not np.isclose(process["shaft_power_W"],load["shaft_power_W"],rtol=1e-9):
                raise RuntimeError("Actual DWSIM pump differs from the Julia layered load.")
            case={"support":support["name"],"julia":load,"dwsim":process,"checks":checks,
                  "actual_specific_shaft_energy_kWh_m3":process["shaft_power_W"]/load["flow_m3_s"]/3.6e6}
            write_json(target / "coupled_result.json",case)
            cases.append(case)
            rows.append({"scenario":name,"alpha_selective_to_baseline_support":alpha,
                "support":support["name"],"radius_m":support["radius_m"],
                "permeability_m2":support["permeability_m2"],"flux_m_s":SETTINGS["flux_m_s"],
                **load,"actual_DWSIM_power_W":process["shaft_power_W"],
                "actual_specific_shaft_energy_kWh_m3":case["actual_specific_shaft_energy_kWh_m3"]})
        saving=1-cases[1]["dwsim"]["shaft_power_W"]/cases[0]["dwsim"]["shaft_power_W"]
        expected=(1-supports[0]["permeability_m2"]/supports[1]["permeability_m2"])/(1+alpha)
        if not np.isclose(saving,expected,rtol=1e-9):
            raise RuntimeError("Layer-weighted saving does not match the independent series-law oracle.")
        scenarios.append({"name":name,"alpha":alpha,"selective_resistance_m_inv":Rl,
            "shaft_energy_saving_fraction":saving,"cases":cases})
    with (directory / "comparison.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    names=("filtration/run.py","filtration/model.py","filtration/build.py","filtration/openlb/capillary3d.cpp",
        "filtration/openlb/config.mk","filtration/openlb/Makefile","hydraulic/openlb/case.h",
        "hydraulic/build.py","hydraulic/run.py","julia/src/Filtration.jl","julia/src/Pore2ProcessDemo.jl",
        "julia/src/Hydraulics.jl","compatibility/run_dwsim.py","compatibility/dwsim_runner/Program.cs",
        "src/pore2process_demo/bridge.py","filtration/check.py","filtration/plot.py","pyproject.toml",
        "uv.lock","julia/Project.toml","julia/Manifest.toml")
    summary={"schema_version":1,"status":"passed","run_id":run_id,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "settings":SETTINGS,"fluid_anchor":anchor,"supports":supports,"validation":extra,"scenarios":scenarios,
        "support_resistance_reduction_fraction":1-supports[0]["permeability_m2"]/supports[1]["permeability_m2"],
        "julia_initialization_s":cold,"reused_pore_run":str(reuse) if reuse else None,
        "experimental_calibration":False,
        "scope":"Synthetic straight-pore support screening, clean water, fixed flux and selective resistance within each pair. No grading, fouling, rejection, compaction or mechanical-strength validation. Total here covers modelled transmembrane hydraulic duty only; shaft energy, without motor, crossflow, backwash or plant auxiliaries. No native DWSIM membrane unit.",
        "units":{"permeability":"m^2","resistance":"m^-1","flux":"m/s","pressure":"Pa","power":"W","specific_shaft_energy":"kWh/m^3"},
        "contracts":"SI scalars; OpenLB +x flow and gauge pressure; 20um square representative cell; original half-cell coordinate/index contract retained. Julia Δp = μJ(Rs+Rl+Rfixed); DWSIM absolute discharge = 1e5 Pa + Δp, mass flow = rho*J*A; eta=0.75 shaft-to-hydraulic.",
        "versions":{"python":sys.version.split()[0],"julia":str(jl.seval("VERSION")),"openlb":"1.9r0","dwsim":anchor["version"],**adapter},
        "vendor_integrity":verify_vendor(),"binary_sha256":sha(BINARY),
        "openlb_build":json.loads((ROOT / "results/filtration/build.json").read_text()),
        "implementation_sha256":{name:sha(ROOT / name) for name in names},
        "runtime_sha256":{str(p.relative_to(ROOT)):sha(p) for p in (ARTIFACTS / "bin/DwsimProbe/debug/DwsimProbe.dll",DISTRIBUTION / "DWSIM.UnitOperations.dll",DISTRIBUTION / "DWSIM.Thermodynamics.dll")}}
    summary["artifacts_sha256"]={str(p.relative_to(directory)):sha(p) for p in directory.rglob("*")
        if p.is_file() and p.suffix in (".csv",".json",".dwxmz")}
    write_json(directory / "summary.json",summary)
    from filtration.plot import render
    render(directory / "summary.json")
    write_json(output / "latest.json",{"run_id":run_id,"summary":str(directory / "summary.json"),"status":"passed"})
    print(f"Fixed-flux filtration study passed: {directory / 'summary.json'}",flush=True)
    return directory / "summary.json"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT / "results/filtration")
    parser.add_argument("--reuse",type=Path,help="Verify and reuse only pore exports from a previous filtration study; rerun DWSIM.")
    args=parser.parse_args()
    build_support()
    execute(args.output.resolve(),args.reuse.resolve() if args.reuse else None)


if __name__=="__main__":
    main()
