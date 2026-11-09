"""Verify saved evidence, then recompute OpenLB analysis and layer loads in Julia.

This reads stored DWSIM results; it does not perform a new DWSIM solve.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from coupling.run import sha, validate_process
from hydraulic.run import analyze, write_json, ROOT
from filtration.model import layer_load


def same(actual,expected,where):
    if actual.keys()!=expected.keys():
        raise RuntimeError(f"Contract fields differ: {where}")
    for name,value in expected.items():
        if isinstance(value,bool):
            ok=actual[name] is value
        else:
            ok=math.isclose(actual[name],value,rel_tol=1e-12,abs_tol=1e-30)
        if not ok:
            raise RuntimeError(f"Recomputed value differs: {where}/{name}")


def verify(path):
    path=Path(path).resolve(); data=json.loads(path.read_text()); base=path.parent
    if data["status"]!="passed" or data["experimental_calibration"] is not False:
        raise RuntimeError("Expected a passed, explicitly synthetic filtration study.")
    for name,digest in data["artifacts_sha256"].items():
        target=(base / name).resolve()
        if base not in target.parents or sha(target)!=digest:
            raise RuntimeError(f"Recorded artifact changed: {name}")
    from pore2process_demo.bridge import runtime
    jl=runtime()
    mu=data["fluid_anchor"]["inlet_viscosity_Pa_s"]; rho=data["fluid_anchor"]["inlet_density_kg_m3"]
    for support in data["supports"]:
        for case in support["pore_cases"]:
            n=case["inputs_and_measurements"]["resolution"]
            directory=base / "supports" / support["name"] / f"n{n}"
            raw=json.loads((directory / "raw.json").read_text())
            if raw!=case["inputs_and_measurements"] or not all(case["checks"].values()):
                raise RuntimeError("OpenLB input snapshot or acceptance differs.")
            same(analyze(directory,raw,jl),case["metrics"],str(directory))
        if not math.isclose(support["permeability_m2"],support["pore_cases"][-1]["metrics"]["permeability_cell_m2"],rel_tol=1e-12):
            raise RuntimeError("Support property differs from its fine-grid measured k.")
    for case in data["validation"]:
        directory=base / "validation" / case["name"]
        raw=json.loads((directory / "raw.json").read_text())
        if raw!=case["inputs_and_measurements"]:
            raise RuntimeError("Validation input snapshot differs.")
        same(analyze(directory,raw,jl),case["metrics"],str(directory))
    settings=data["settings"]
    if not math.isclose(settings["flux_m_s"],settings["flux_LMH"]*1e-3/3600,rel_tol=1e-14):
        raise RuntimeError("The declared filtration flux violates the LMH → SI contract.")
    with (base / "comparison.csv").open(newline="") as stream:
        csv_rows=list(csv.DictReader(stream))
    expected_rows=sum(len(s["cases"]) for s in data["scenarios"])
    if len(csv_rows)!=expected_rows:
        raise RuntimeError("Comparison CSV has an unexpected row count.")
    row=0
    Rs0=settings["support_thickness_m"]/data["supports"][0]["permeability_m2"]
    for scenario in data["scenarios"]:
        if not math.isclose(scenario["selective_resistance_m_inv"],scenario["alpha"]*Rs0,rel_tol=1e-12):
            raise RuntimeError("Selective resistance differs from its scenario definition.")
        for case,support in zip(scenario["cases"],data["supports"],strict=True):
            directory=base / "scenarios" / scenario["name"] / case["support"]
            process=json.loads((directory / "summary.json").read_text())
            inputs=json.loads((directory / "inputs.json").read_text())
            if process!=case["dwsim"] or process["inputs"]!=inputs or case["support"]!=support["name"]:
                raise RuntimeError("DWSIM, input or support snapshots differ.")
            if case!=json.loads((directory / "coupled_result.json").read_text()):
                raise RuntimeError("Coupled result file differs from the main snapshot.")
            expected=layer_load(support["permeability_m2"],settings,scenario["selective_resistance_m_inv"],mu)
            same(case["julia"],expected,"layer load")
            if not math.isclose(inputs["Outlet_pressure_Pa"],settings["inlet_absolute_pressure_Pa"]+expected["total_pressure_drop_Pa"],rel_tol=1e-12):
                raise RuntimeError("DWSIM discharge does not contain the Julia layered load.")
            if not math.isclose(inputs["Mass_flow_kg_s"],rho*expected["flow_m3_s"],rel_tol=1e-12):
                raise RuntimeError("DWSIM mass flow differs from fixed flux times membrane area.")
            validate_process(process,inputs,rho,mu,expected["total_pressure_drop_Pa"])
            energy=process["shaft_power_W"]/expected["flow_m3_s"]/3.6e6
            if not math.isclose(energy,case["actual_specific_shaft_energy_kWh_m3"],rel_tol=1e-12):
                raise RuntimeError("Exported specific energy differs from actual pump power.")
            csv_expected={"alpha_selective_to_baseline_support":scenario["alpha"],
                "radius_m":support["radius_m"],"permeability_m2":support["permeability_m2"],
                "flux_m_s":settings["flux_m_s"],**expected,"actual_DWSIM_power_W":process["shaft_power_W"],
                "actual_specific_shaft_energy_kWh_m3":energy}
            if csv_rows[row]["scenario"]!=scenario["name"] or csv_rows[row]["support"]!=support["name"]:
                raise RuntimeError("CSV scenario or support identification differs.")
            for field,value in csv_expected.items():
                exported=csv_rows[row][field]
                if isinstance(value,bool):
                    ok=exported==str(value)
                else:
                    ok=math.isclose(float(exported),value,rel_tol=1e-12,abs_tol=1e-30)
                if not ok:
                    raise RuntimeError(f"CSV differs from Julia: {field}")
            row+=1
        saving=1-scenario["cases"][1]["dwsim"]["shaft_power_W"]/scenario["cases"][0]["dwsim"]["shaft_power_W"]
        independent=(1-data["supports"][0]["permeability_m2"]/data["supports"][1]["permeability_m2"])/(1+scenario["alpha"])
        if not math.isclose(saving,scenario["shaft_energy_saving_fraction"],rel_tol=1e-12) or not math.isclose(saving,independent,rel_tol=1e-9):
            raise RuntimeError("The saved energy saving differs from series-resistance theory.")
    result={"status":"passed","summary":str(path),"verified_artifacts":len(data["artifacts_sha256"]),
            "recomputed_Julia_loads":row,"actual_new_DWSIM_solves":0}
    write_json(base / "verification.json",result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("summary",type=Path,nargs="?")
    path=parser.parse_args().summary
    if path is None:
        path=Path(json.loads((ROOT / "results/filtration/latest.json").read_text())["summary"])
    print(json.dumps(verify(path),indent=2))
