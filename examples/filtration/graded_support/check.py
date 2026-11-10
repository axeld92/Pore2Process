"""Replay stored native analysis and Julia loads; inspect saved DWSIM evidence.

This does not rerun OpenLB or DWSIM. Their actual solve evidence is hash-checked.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from examples.filtration.calibration.check import close_tree
from pore2process_demo.io import sha; from pore2process_demo.runtime.hydraulics import validate_process
from pore2process_demo.io import write_json
from examples.filtration.graded_support.geometry import audit
from examples.filtration.graded_support.model import ROOT, analyze, load_case

def verify(path, *, write_record=True):
    path=Path(path).resolve(); base=path.parent
    s=json.loads(path.read_text()); p=s["protocol"]
    if s["status"]!="passed": raise RuntimeError("Expected a completed study.")
    actual_files={str(f.relative_to(base)) for f in base.rglob("*") if f.is_file() and str(f.relative_to(base)) not in ("summary.json","status.json","verification.json")}
    if set(s["artifact_sha256"])!=actual_files:
        raise RuntimeError("Artifact coverage differs from saved directory.")
    for name,digest in s["artifact_sha256"].items():
        f=(base/name).resolve()
        if base not in f.parents or sha(f)!=digest:
            raise RuntimeError(f"Artifact hash mismatch: {name}")
    inputs=json.loads((base/"inputs.json").read_text())
    for key in ("protocol","literature_reference","fluid_anchor","binary_sha256","implementation_sha256"):
        close_tree(s[key],inputs[key],key)
    from pore2process_demo.bridge import runtime
    jl=runtime()
    mu=s["fluid_anchor"]["inlet_viscosity_Pa_s"]; rho=s["fluid_anchor"]["inlet_density_kg_m3"]
    from examples.filtration.graded_support.run import requests, MEASUREMENTS
    jobs,names=requests(p,mu,rho)
    if set(s["native_cases"])!={name for name,_ in jobs}:
        raise RuntimeError("Missing/extra native validation cases.")
    for name,requested in jobs:
        record=s["native_cases"][name]; directory=base/name
        measurement=json.loads((directory/"measurement.json").read_text())
        close_tree(record["measurement"],measurement,"measurement")
        close_tree(measurement["inputs"],requested,"requested native inputs")
        if measurement["binary_sha256"]!=s["binary_sha256"] or set(measurement["artifact_sha256"])!=set(MEASUREMENTS):
            raise RuntimeError("Native binary/file contract mismatch.")
        for filename,digest in measurement["artifact_sha256"].items():
            if sha(directory/filename)!=digest: raise RuntimeError("Native measurement changed.")
        raw=json.loads((directory/"raw.json").read_text())
        close_tree(record["raw"],raw,"raw")
        for field,value in (("resolution",requested["resolution"]),("bulk_radius_m",requested["bulk_radius_m"]),
                ("lead_length_each_m",requested["lead_m"]),("tau",requested["tau"]),
                ("imposed_delta_p_Pa",requested["pressure_Pa"]),("viscosity_Pa_s",mu),("reference_density_kg_m3",rho),
                ("body_length_m",p["geometry"]["body_length_m"]),
                ("cell_width_m",p["geometry"]["representative_cell_width_m"])):
            close_tree(raw[field],value,f"native.{field}")
        metrics,checks=analyze(directory,raw,jl,p)
        close_tree(record["metrics"],metrics,"native analysis")
        close_tree(record["checks"],checks,"native checks")
    if [v["name"] for v in s["supports"]]!=list(names):
        raise RuntimeError("Geometry ordering differs.")
    for support in s["supports"]:
        name=support["name"]; grids=support["grids"]
        for grid,n in zip(grids,p["numerics"]["resolutions_across_neck_diameter"],strict=True):
            close_tree(grid,s["native_cases"][f"supports/{name}/n{n}"],"grid")
        R=grids[-1]["metrics"]["body_resistance_Pa_s_m3"]
        close_tree(support["bulk_radius_m"],grids[-1]["raw"]["bulk_radius_m"],"bulk radius")
        nominal_flow=p["process"]["filtration_flux_LMH"]/3.6e6*p["geometry"]["representative_cell_width_m"]**2
        nominal_Re=grids[-1]["metrics"]["max_Re"]*nominal_flow/abs(grids[-1]["metrics"]["measured_flow_m3_s"])
        close_tree(support["uncalibrated_geometry_nominal_flux_Re"],nominal_Re,"uncalibrated geometry nominal-flux Re")
        if nominal_Re>=p["acceptance"]["native_max_Re"]: raise RuntimeError("Nominal-flux Re screen failed.")
        change=abs(R/grids[-2]["metrics"]["body_resistance_Pa_s_m3"]-1)
        close_tree(support["body_resistance_Pa_s_m3"],R,"measured R")
        close_tree(support["fine_grid_relative_change"],change,"grid change")
        if change>=p["acceptance"]["fine_grid_body_resistance_relative_change"] or not all(support["checks"].values()):
            raise RuntimeError("Grid/cross-case acceptance failed.")
        expected_ids={case_id for case_id,_ in jobs if case_id.startswith(f"validation/{name}/")}
        if {v["case_id"] for v in support["validation"]}!=expected_ids or len(support["validation"])!=len(expected_ids):
            raise RuntimeError("Cross-case validation coverage differs.")
        expected_checks={"fine_grid_resistance":change<p["acceptance"]["fine_grid_body_resistance_relative_change"]}
        for v in support["validation"]:
            value=abs(s["native_cases"][v["case_id"]]["metrics"]["body_resistance_Pa_s_m3"]/grids[0]["metrics"]["body_resistance_Pa_s_m3"]-1)
            close_tree(value,v["relative_resistance_change"],"cross-case change")
            key="lead_extension_resistance_relative_change" if v["case_id"].endswith("extended_leads") else "tau_resistance_relative_change" if v["case_id"].endswith("tau_check") else "pressure_range_resistance_relative_change"
            close_tree(v["criterion"],p["acceptance"][key],"frozen criterion")
            if value>=p["acceptance"][key]: raise RuntimeError("Cross-case screen failed.")
            expected_checks[v["case_id"]]=value<p["acceptance"][key]
        close_tree(support["checks"],expected_checks,"cross-case acceptance flags")
        close_tree(audit(base/grids[-1]["case_id"]/"geometry/voxels.csv",grids[-1]["raw"],p),support["geometry_connectivity"],"connectivity")
    R0=s["supports"][0]["body_resistance_Pa_s_m3"]
    source=s["literature_reference"]
    lo,hi=source["digitization_sensitivity"]["fixed_impedance_range_Pa_s_m"]
    scenarios=(("nominal",source["fit"]["fixed_impedance_Pa_s_m"]),("reading_low",lo),("reading_high",hi))
    if [v["name"] for v in s["scenarios"]]!=[name for name,_ in scenarios]:
        raise RuntimeError("Reading scenarios differ.")
    with (base/"comparison.csv").open(newline="") as f:
        csv_rows=list(csv.DictReader(f))
    if len(csv_rows)!=9: raise RuntimeError("Expected nine process rows.")
    counter=0
    for scenario,(label,Z) in zip(s["scenarios"],scenarios,strict=True):
        close_tree(scenario["fixed_non_support_impedance_Pa_s_m"],Z,"fixed non-support Z")
        cases=scenario["cases"]
        for case,support in zip(cases,s["supports"],strict=True):
            close_tree(case["support"],support["name"],"process support ID")
            expected=load_case(jl,support["body_resistance_Pa_s_m3"],R0,p,Z)
            close_tree(case["julia"],expected,"Julia scale translation")
            directory=base/"process"/label/support["name"]
            actual=json.loads((directory/"summary.json").read_text())
            actual_inputs=json.loads((directory/"inputs.json").read_text())
            close_tree(actual,case["dwsim"],"DWSIM saved evidence")
            close_tree(actual["inputs"],actual_inputs,"DWSIM input snapshot")
            for field in ("Temperature_K","Pump_efficiency_fraction"):
                key="temperature_K" if field=="Temperature_K" else "pump_efficiency_fraction"
                close_tree(actual_inputs[field],p["process"][key],f"DWSIM fixed {field}")
            close_tree(actual_inputs["Inlet_pressure_Pa"],p["process"]["inlet_absolute_pressure_Pa"],"DWSIM fixed inlet pressure")
            close_tree(actual_inputs["Outlet_pressure_Pa"],p["process"]["inlet_absolute_pressure_Pa"]+expected["pressure_drop_Pa"],"pressure supplied to DWSIM")
            close_tree(actual_inputs["Mass_flow_kg_s"],rho*expected["flow_m3_s"],"mass flow supplied to DWSIM")
            checks=validate_process(actual,actual_inputs,rho,mu,expected["pressure_drop_Pa"])
            close_tree(checks,case["checks"],"DWSIM acceptance")
            if not (directory/"water_pump.dwxmz").is_file(): raise RuntimeError("Missing actual saved flowsheet.")
            energy=actual["shaft_power_W"]/expected["flow_m3_s"]/3.6e6
            close_tree(case["actual_specific_shaft_energy_kWh_m3"],energy,"specific energy")
            saving=1-actual["shaft_power_W"]/cases[0]["dwsim"]["shaft_power_W"]
            independent=cases[0]["julia"]["support_fraction"]*(1-expected["measured_body_resistance_ratio"])
            close_tree(case["energy_saving_vs_straight_fraction"],saving,"actual design saving")
            close_tree(saving,independent,"independent layer-weighted energy oracle")
            row=csv_rows[counter]; counter+=1
            for field,value in {"scenario":label,"support":support["name"],"bulk_radius_m":support["bulk_radius_m"],
                    **expected,"actual_DWSIM_power_W":actual["shaft_power_W"],"actual_specific_shaft_energy_kWh_m3":energy}.items():
                close_tree(row[field] if isinstance(value,str) else float(row[field]),value,f"CSV.{field}")
    verification={"status":"passed","native_cases_reanalyzed":len(jobs),"connectivity_audits":3,
        "actual_saved_DWSIM_cases_checked":9,"Julia_loads_replayed":9,"artifacts_checked":len(s["artifact_sha256"]),
        "new_OpenLB_or_DWSIM_solves_performed_by_checker":False}
    if write_record:
        write_json(base/"verification.json",verification)
    return verification

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument("summary",nargs="?",type=Path)
    args=a.parse_args()
    path=args.summary or Path(json.loads((ROOT/"results/graded_support/latest.json").read_text())["summary"])
    print(verify(path))

if __name__=="__main__":
    main()
