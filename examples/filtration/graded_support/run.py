"""Three actual OpenLB geometries → measured ratios in Julia → DWSIM pump cases."""
from pore2process_demo.paths import project_path
from pore2process_demo.provenance import implementation_files, source_manifest
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np
from pore2process_demo.runtime.dwsim import build_adapter, solve_case, ARTIFACTS, DISTRIBUTION
from pore2process_demo.io import sha; from pore2process_demo.runtime.hydraulics import validate_process
from pore2process_demo.runtime.openlb import verify_vendor
from pore2process_demo.runtime.hydraulics import simulate; from pore2process_demo.io import write_json
from examples.filtration.graded_support.build import BINARY, main as build
from examples.filtration.graded_support.geometry import audit
from examples.filtration.graded_support.model import ROOT, read_protocol, analyze, load_case

MEASUREMENTS=("raw.json","planes.csv","midplane.csv","history.csv","stdout.log")

def native_case(base,case_id,inputs,protocol,reuse=None):
    """Workers only launch isolated native processes; Julia stays on main thread."""
    target=base/case_id
    signature={"binary_sha256":sha(BINARY),"inputs":inputs}
    if reuse and (reuse/case_id/"measurement.json").exists():
        origin=reuse/case_id
        record=json.loads((origin/"measurement.json").read_text())
        if {k:record[k] for k in signature}!=signature or set(record["artifact_sha256"])!=set(MEASUREMENTS):
            raise RuntimeError("Native cache binary/input/artifact contract differs.")
        target.mkdir(parents=True)
        for name,digest in record["artifact_sha256"].items():
            if sha(origin/name)!=digest:
                raise RuntimeError("Native cache artifact changed.")
            shutil.copy2(origin/name,target/name)
        raw=json.loads((target/"raw.json").read_text())
        command=record["command"]
        reused=str(origin)
    else:
        n=protocol["numerics"]
        raw,command=simulate(target,inputs["resolution"],pressure=inputs["pressure_Pa"],
            tau=inputs["tau"],tolerance=n["steady_tolerance"],radius=protocol["geometry"]["neck_radius_m"],
            viscosity=inputs["viscosity_Pa_s"],density=inputs["density_kg_m3"],binary=BINARY,
            timeout_s=3600,
            options=("--bulk-radius-m",str(inputs["bulk_radius_m"]),"--lead-length-m",str(inputs["lead_m"])))
        reused=None
    for key,value in (("bulk_radius_m",inputs["bulk_radius_m"]),("lead_length_each_m",inputs["lead_m"]),
                      ("body_length_m",protocol["geometry"]["body_length_m"]),
                      ("cell_width_m",protocol["geometry"]["representative_cell_width_m"])):
        if raw[key]!=value:
            raise RuntimeError(f"Native geometry snapshot differs: {key}")
    record={**signature,"case_id":case_id,"command":command,"reused_native":reused,
            "artifact_sha256":{name:sha(target/name) for name in MEASUREMENTS}}
    write_json(target/"measurement.json",record)
    return record,raw

def requests(protocol,mu,rho):
    g,n=protocol["geometry"],protocol["numerics"]
    jobs=[]
    names=("straight","widened_5um","widened_6um")
    for name,bulk in zip(names,g["bulk_radii_m"],strict=True):
        common={"bulk_radius_m":bulk,"lead_m":g["lead_length_each_m"],
                "tau":n["tau"],"viscosity_Pa_s":mu,"density_kg_m3":rho,"pressure_Pa":n["pressure_drop_Pa"]}
        for grid in n["resolutions_across_neck_diameter"]:
            jobs.append((f"supports/{name}/n{grid}",{**common,"resolution":grid}))
        for label,overrides in (
            *((f"pressure_{dp:g}",{"pressure_Pa":dp}) for dp in n["pressure_checks_Pa"]),
            ("extended_leads",{"lead_m":g["domain_check_lead_length_each_m"]}),
            ("tau_check",{"tau":n["tau_check"]})):
            jobs.append((f"validation/{name}/{label}",{**common,"resolution":n["resolutions_across_neck_diameter"][0],**overrides}))
    return jobs,names

def execute(output,reuse=None,workers=2,protocol_path=None):
    p=read_protocol(protocol_path); g,n,limits=p["geometry"],p["numerics"],p["acceptance"]
    reference_path=project_path(ROOT/p["process"]["literature_reference"])
    reference=json.loads(reference_path.read_text())
    if reference["protocol_id"]!="zhu-2021-M1-series-impedance-pdf-v2":
        raise RuntimeError("Review a changed literature reference before running.")
    run_id=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base=output/"runs"/run_id; base.mkdir(parents=True)
    adapter=build_adapter(base/"adapter_build")
    process=p["process"]
    anchor_inputs={"Temperature_K":process["temperature_K"],"Inlet_pressure_Pa":process["inlet_absolute_pressure_Pa"],
        "Outlet_pressure_Pa":process["inlet_absolute_pressure_Pa"]+1000.,"Mass_flow_kg_s":.01,
        "Pump_efficiency_fraction":process["pump_efficiency_fraction"],"Scope":"DWSIM liquid-water property anchor for graded-support scenarios."}
    anchor=solve_case(base/"fluid_anchor",anchor_inputs)
    mu,rho=anchor["inlet_viscosity_Pa_s"],anchor["inlet_density_kg_m3"]
    if anchor["inlet_vapor_fraction"]!=0 or min(mu,rho)<=0:
        raise RuntimeError("Shared property anchor must be liquid water with positive viscosity/density.")
    from pore2process_demo.bridge import runtime,runtime_versions
    started=perf_counter(); jl=runtime(); cold=perf_counter()-started
    code_paths=sorted([f for f in (ROOT/"examples/filtration/graded_support").rglob("*") if f.is_file() and 'reference' not in f.parts and f.suffix in (".py",".h",".cpp",".mk",".json")]+[ROOT/"examples/filtration/graded_support/openlb/Makefile"])
    code_paths += [project_path(ROOT/f) for f in ("julia/src/GradedSupport.jl","julia/src/Pore2ProcessDemo.jl",
        "julia/src/SeriesCalibration.jl","julia/src/Hydraulics.jl","examples/hydraulics/capillary_benchmark/build.py","examples/hydraulics/capillary_benchmark/run.py",
        "examples/hydraulics/darcy_pump/run.py","src/pore2process_demo/runtime/dwsim.py","src/pore2process_demo/runtime/dwsim_runner/Program.cs",
        "src/pore2process_demo/bridge.py","uv.lock","julia/Project.toml","julia/Manifest.toml")]
    code_paths = implementation_files(code_paths)
    snapshot={"protocol":p,"literature_reference":reference,"literature_reference_sha256":sha(reference_path),
        "fluid_anchor":anchor,"binary_sha256":sha(BINARY),
        "implementation_sha256":source_manifest(code_paths)}
    write_json(base/"inputs.json",snapshot)
    write_json(base/"status.json",{"status":"running","run_id":run_id})
    jobs,names=requests(p,mu,rho)
    records={}
    pool=ThreadPoolExecutor(max_workers=workers)
    try:
        futures={pool.submit(native_case,base,name,inputs,p,reuse):name for name,inputs in jobs}
        for future in as_completed(futures):
            name=futures[future]
            measurement,raw=future.result()
            metrics,checks=analyze(base/name,raw,jl,p)
            record={"case_id":name,"raw":raw,"metrics":metrics,"checks":checks,"measurement":measurement}
            records[name]=record
            write_json(base/name/"analysis.json",record)
            print(f"Accepted {name}: Rbody={metrics['body_resistance_Pa_s_m3']:.8g}; mass spread={metrics['all_plane_mass_relative_spread']:.3%}",flush=True)
    except BaseException as exc:
        write_json(base/"status.json",{"status":"failed","reason":str(exc),"completed_analysis":list(records)})
        pool.shutdown(wait=True,cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    supports=[]
    for name in names:
        grids=[records[f"supports/{name}/n{grid}"] for grid in n["resolutions_across_neck_diameter"]]
        fine=grids[-1]; R=fine["metrics"]["body_resistance_Pa_s_m3"]
        grid_change=abs(R/grids[-2]["metrics"]["body_resistance_Pa_s_m3"]-1)
        checks={"fine_grid_resistance":grid_change<limits["fine_grid_body_resistance_relative_change"]}
        validation=[]
        for case_id,record in records.items():
            if not case_id.startswith(f"validation/{name}/"): continue
            change=abs(record["metrics"]["body_resistance_Pa_s_m3"]/grids[0]["metrics"]["body_resistance_Pa_s_m3"]-1)
            key="lead_extension_resistance_relative_change" if case_id.endswith("extended_leads") else "tau_resistance_relative_change" if case_id.endswith("tau_check") else "pressure_range_resistance_relative_change"
            checks[case_id]=change<limits[key]
            validation.append({"case_id":case_id,"relative_resistance_change":change,"criterion":limits[key]})
        if not all(checks.values()):
            write_json(base/"status.json",{"status":"failed","reason":"cross-case checks","support":name,"checks":checks,"validation":validation,"grid_relative_change":grid_change})
            raise RuntimeError(f"Cross-case acceptance failed for {name}: {checks}; inspect {base}")
        geometry=base/fine["case_id"] / "geometry"; geometry.mkdir()
        command=[str(BINARY),"--geometry-only","--resolution",str(fine["raw"]["resolution"]),
                 "--bulk-radius-m",str(fine["raw"]["bulk_radius_m"])]
        with (geometry/"stdout.log").open("w") as log:
            subprocess.run(command,cwd=geometry,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
        connectivity=audit(geometry/"voxels.csv",fine["raw"],p)
        nominal_flow=process["filtration_flux_LMH"]/3.6e6*g["representative_cell_width_m"]**2
        flux_scale=nominal_flow/abs(fine["metrics"]["measured_flow_m3_s"])
        nominal_Re=fine["metrics"]["max_Re"]*flux_scale
        if nominal_Re>=limits["native_max_Re"]:
            raise RuntimeError("Uncalibrated geometry at nominal flux leaves the checked creeping-flow range.")
        supports.append({"name":name,"bulk_radius_m":fine["raw"]["bulk_radius_m"],"grids":grids,
            "body_resistance_Pa_s_m3":R,"fine_grid_relative_change":grid_change,"checks":checks,
            "validation":validation,"geometry_connectivity":connectivity,
            "uncalibrated_geometry_nominal_flux_Re":nominal_Re})
    R0=supports[0]["body_resistance_Pa_s_m3"]
    if not R0>supports[1]["body_resistance_Pa_s_m3"]>supports[2]["body_resistance_Pa_s_m3"]:
        raise RuntimeError("Widening did not produce the expected resolved resistance ordering.")
    Z=reference["fit"]["fixed_impedance_Pa_s_m"]
    lo,hi=reference["digitization_sensitivity"]["fixed_impedance_range_Pa_s_m"]
    scenarios=[]; rows=[]
    for label,impedance in (("nominal",Z),("reading_low",lo),("reading_high",hi)):
        cases=[]
        for support in supports:
            load=load_case(jl,support["body_resistance_Pa_s_m3"],R0,p,impedance)
            inputs={"Temperature_K":process["temperature_K"],"Inlet_pressure_Pa":process["inlet_absolute_pressure_Pa"],
                "Outlet_pressure_Pa":process["inlet_absolute_pressure_Pa"]+load["pressure_drop_Pa"],
                "Mass_flow_kg_s":rho*load["flow_m3_s"],"Pump_efficiency_fraction":process["pump_efficiency_fraction"],
                "Scope":f"Graded-support {run_id}/{label}/{support['name']}. Measured OpenLB resistance ratio; common literature baseline calibration; fixed effective non-support Z. Transmembrane load only, no native membrane separation unit."}
            directory=base/"process"/label/support["name"]
            actual=solve_case(directory,inputs)
            checks=validate_process(actual,inputs,rho,mu,load["pressure_drop_Pa"])
            if not np.isclose(actual["shaft_power_W"],load["shaft_power_W"],rtol=limits["dwsim_power_relative_error"]):
                raise RuntimeError("DWSIM power differs from the Julia hydraulic load.")
            case={"support":support["name"],"julia":load,"dwsim":actual,"checks":checks,
                "actual_specific_shaft_energy_kWh_m3":actual["shaft_power_W"]/load["flow_m3_s"]/3.6e6}
            cases.append(case)
            rows.append({"scenario":label,"support":support["name"],"bulk_radius_m":support["bulk_radius_m"],
                **load,"actual_DWSIM_power_W":actual["shaft_power_W"],
                "actual_specific_shaft_energy_kWh_m3":case["actual_specific_shaft_energy_kWh_m3"]})
            write_json(directory/"coupled_result.json",case)
        for case in cases:
            case["energy_saving_vs_straight_fraction"]=1-case["dwsim"]["shaft_power_W"]/cases[0]["dwsim"]["shaft_power_W"]
        scenarios.append({"name":label,"fixed_non_support_impedance_Pa_s_m":impedance,"cases":cases})
    import csv
    with (base/"comparison.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary={"schema_version":1,"status":"passed","run_id":run_id,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "protocol":p,"fluid_anchor":anchor,"supports":supports,"native_cases":records,"scenarios":scenarios,
        "literature_reference":reference,"literature_reference_sha256":sha(reference_path),
        "julia_initialization_s":cold,"versions":runtime_versions(),"adapter":adapter,"native_workers":workers,
        "binary_sha256":sha(BINARY),"implementation_sha256":snapshot["implementation_sha256"],
        "runtime_sha256":{str(f.relative_to(ROOT)):sha(f) for f in (ARTIFACTS/"bin/DwsimProbe/debug/DwsimProbe.dll",DISTRIBUTION/"DWSIM.UnitOperations.dll",DISTRIBUTION/"DWSIM.Thermodynamics.dll")},
        "vendor_integrity":verify_vendor(),"build":json.loads((ROOT/"results/graded_support/build.json").read_text()),
        "units":{"length":"m","pressure":"Pa gauge for OpenLB; absolute for DWSIM","flow":"m3/s, positive +x","native_body_resistance":"Pa s/m3","impedance":"Pa s/m","permeability":"m2","shaft_power":"W","specific_shaft_energy":"kWh/m3"},
        "scope":process["scope"],"translation_assumption":process["translation"]}
    write_json(base/"summary.json",summary)
    from examples.filtration.graded_support.plot import render
    render(base/"summary.json")
    summary["artifact_sha256"]={str(f.relative_to(base)):sha(f) for f in sorted(base.rglob("*")) if f.is_file() and str(f.relative_to(base)) not in ("summary.json","status.json","verification.json")}
    write_json(base/"summary.json",summary)
    from examples.filtration.graded_support.check import verify
    try:
        verify(base/"summary.json")
    except BaseException as exc:
        summary["status"]="failed_verification"
        write_json(base/"summary.json",summary)
        write_json(base/"status.json",{"status":"failed","reason":str(exc)})
        raise
    write_json(base/"status.json",{"status":"passed","run_id":run_id})
    write_json(output/"latest.json",{"status":"passed","summary":str(base/"summary.json"),"run_id":run_id})
    print(f"Graded-support study passed: {base/'summary.json'}",flush=True)
    return base/"summary.json"

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/graded_support")
    parser.add_argument("--reuse-native",type=Path,help="Reuse only hash-verified native measurements from a prior run directory, with identical binary/inputs; rerun analysis and DWSIM.")
    parser.add_argument("--workers",type=int,choices=(1,2),default=2)
    parser.add_argument("--protocol",type=Path,help="Explicit versioned protocol; saved in each run.")
    args=parser.parse_args()
    build()
    execute(args.output.resolve(),args.reuse_native.resolve() if args.reuse_native else None,args.workers,args.protocol)

if __name__=="__main__":
    main()
