"""Audited literature inputs → real Julia log fit → held-out diagnostics/exports.

Run from the repository root: uv run python -m calibration.run
No geometry-resolved or complete-membrane experimental validation is implied.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter
from uuid import uuid4

from pore2process_demo.bridge import ROOT, runtime, runtime_versions
from filtration.model import layer_load

ROLE_KEYS = {"calibration": "calibration_rows", "primary_holdout": "primary_holdout_rows",
             "transfer": "transfer_rows", "out_of_domain": "out_of_domain_rows"}
METRICS = ("count", "mean_absolute_relative_error", "max_absolute_relative_error",
           "mean_relative_error", "log_rmse")


def read_json(path):
    def reject(value):
        raise ValueError(f"Nonstandard JSON value: {value}")
    return json.loads(Path(path).read_text(), parse_constant=reject)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def roles_for(rows, protocol):
    ids = [r["table_row"] for r in rows]
    assigned = [n for key in ROLE_KEYS.values() for n in protocol[key]]
    if len(set(ids)) != len(ids) or len(set(assigned)) != len(assigned) or set(ids) != set(assigned):
        raise ValueError("Protocol split must cover each unique source row exactly once, without overlap.")
    return {n: role for role, key in ROLE_KEYS.items() for n in protocol[key]}


def validate_snapshot(snapshot):
    p, audit, source = snapshot["protocol"], snapshot["audit"], snapshot["source"]
    if p["dataset_id"] != audit["dataset_id"]:
        raise ValueError("Protocol and audited dataset differ.")
    rows = source["table2_rows"]
    roles_for(rows,p)
    if p["nominal_primary_particle_diameter_m"] != audit["nominal_primary_particle_diameter_m"]:
        raise ValueError("Protocol particle sizes differ from the source audit.")
    audited = audit["visually_transcribed_table2_rows"]
    if len(rows) != len(audited):
        raise ValueError("Source row count differs from the visual audit.")
    names=("powder", "reported_sinter_temperature_C", "reported_sinter_hold_time_h",
           "reported_k_mantissa_times_1e_minus16_m2", "porosity_fraction", "reported_mean_strength_MPa", "weibull_modulus")
    for expected_id,(row, reference) in enumerate(zip(rows,audited,strict=True),start=1):
        if row["table_row"] != expected_id:
            raise ValueError("Source row IDs must retain Table 2 order.")
        if [row[k] for k in names] != reference:
            raise ValueError(f"Row {row['table_row']} differs from the visually audited table.")
        if not math.isclose(row["permeability_m2"],reference[3]*1e-16,rel_tol=1e-14):
            raise ValueError("Permeability SI multiplier does not match Table 2.")
    limit = p["primary_screen"]["max_absolute_relative_error"]
    if isinstance(limit,bool) or not math.isfinite(limit) or limit <= 0:
        raise ValueError("Primary screening tolerance must be finite and positive.")


def load_snapshot(require_sources=False):
    protocol=read_json(ROOT/"calibration/protocol.json")
    extraction=read_json(ROOT/protocol["extraction_path"])
    source=extraction["biesheuvel_1999"]
    audit=read_json(ROOT/protocol["audit_path"])
    statuses=[]
    hashes={str(Path(source["source"]["folder"])/name):digest
            for name,digest in source["source"]["files_sha256"].items()}
    hashes.update(audit["files_sha256"])
    for name,expected in hashes.items():
        path=ROOT/name
        if not path.exists():
            if require_sources:
                raise ValueError(f"Required local source is missing: {name}")
            status="missing; using retained audited extraction"
        else:
            if sha256(path) != expected:
                raise ValueError(f"Source hash changed; review before fitting: {name}")
            status="sha256 verified"
        statuses.append({"path":name,"sha256":expected,"status":status})
    tracked=["calibration/protocol.json",protocol["extraction_path"],protocol["audit_path"],
             "calibration/run.py","calibration/plot.py","calibration/check.py",
             "filtration/model.py","julia/src/Calibration.jl","julia/src/Filtration.jl",
             "julia/src/Pore2ProcessDemo.jl","src/pore2process_demo/bridge.py",
             "pyproject.toml","uv.lock","julia/Project.toml","julia/Manifest.toml"]
    snapshot={"protocol":protocol,"source":source,"audit":audit,"source_files":statuses,
              "implementation_sha256":{name:sha256(ROOT/name) for name in tracked}}
    validate_snapshot(snapshot)
    return snapshot


def analyze(snapshot):
    """Selection in Python; fit, predictions, error metrics and loads in Julia."""
    validate_snapshot(snapshot)
    p=snapshot["protocol"]
    rows=snapshot["source"]["table2_rows"]
    roles=roles_for(rows,p)
    started=perf_counter()
    jl=runtime()
    initialization_s=perf_counter()-started
    started=perf_counter()
    vector=jl.seval("Vector{Float64}")
    def fit(selected):
        return jl.fit_support_multiplier(
            vector([p["nominal_primary_particle_diameter_m"][r["powder"]] for r in selected]),
            vector([r["porosity_fraction"] for r in selected]),
            vector([r["permeability_m2"] for r in selected]))
    training=[r for r in rows if roles[r["table_row"]]=="calibration"]
    fitted=fit(training)
    b=float(fitted.multiplier)
    cases=[]
    for row in rows:
        d=p["nominal_primary_particle_diameter_m"][row["powder"]]
        base=float(jl.support_permeability(d,row["porosity_fraction"]))
        pred=float(jl.support_permeability(d,row["porosity_fraction"],multiplier=b))
        measured=row["permeability_m2"]
        cases.append({"table_row":row["table_row"],"powder":row["powder"],"role":roles[row["table_row"]],
            "porosity_fraction":row["porosity_fraction"],"nominal_primary_particle_diameter_m":d,
            "measured_permeability_m2":measured,"baseline_CK_permeability_m2":base,
            "predicted_permeability_m2":pred,"relative_error":pred/measured-1,
            "baseline_relative_error":base/measured-1,"permeability_uncertainty_m2":None})
    metrics={}
    for role in ROLE_KEYS:
        group=[c for c in cases if c["role"]==role]
        metrics[role]={}
        for model,key in (("baseline_CK","baseline_CK_permeability_m2"),("calibrated","predicted_permeability_m2")):
            values=jl.permeability_errors(vector([c[key] for c in group]),
                                         vector([c["measured_permeability_m2"] for c in group]))
            metrics[role][model]={name:int(getattr(values,name)) if name=="count" else float(getattr(values,name)) for name in METRICS}
    loo=[]
    for excluded in training:
        retained=[r for r in training if r["table_row"] != excluded["table_row"]]
        value=fit(retained)
        loo.append({"omitted_calibration_row":excluded["table_row"],"multiplier":float(value.multiplier)})
    limit=p["primary_screen"]["max_absolute_relative_error"]
    failed=[c["table_row"] for c in cases if c["role"]=="primary_holdout" and abs(c["relative_error"])>limit]
    process=[]
    settings=p["conditional_process_demo"]
    by_id={c["table_row"]:c for c in cases}
    base,alt=by_id[settings["baseline_row"]],by_id[settings["alternative_row"]]
    reference_R=settings["support_thickness_m"]/base["measured_permeability_m2"]
    for ratio in settings["selective_to_baseline_support_resistance_ratios"]:
        # The same absolute Rselective in measured and predicted calculations.
        Rselective=ratio*reference_R
        cell_max=snapshot["source"]["measurement_conditions"]["reported_cell_pressurization_up_to_Pa"]
        record={"selective_to_measured_baseline_support_ratio":ratio,"fixed_selective_resistance_m_inv":Rselective,
                "reported_experimental_cell_pressurization_max_Pa":cell_max,
                "pressure_limit_scope":"Reported apparatus pressurization only, not specimen-specific Darcy validity or allowable device pressure."}
        for model,key in (("measured_support_k","measured_permeability_m2"),("calibrated_support_k","predicted_permeability_m2")):
            first=layer_load(base[key],settings,Rselective,settings["viscosity_Pa_s"])
            second=layer_load(alt[key],settings,Rselective,settings["viscosity_Pa_s"])
            record[model]={"baseline":first,"alternative":second,
                "any_total_pressure_exceeds_reported_cell_pressurization":max(first["total_pressure_drop_Pa"],second["total_pressure_drop_Pa"])>cell_max,
                "fractional_shaft_energy_saving":1-second["specific_shaft_energy_kWh_m3"]/first["specific_shaft_energy_kWh_m3"]}
        process.append(record)
    return {"fit":{"multiplier":b,"log_multiplier":float(fitted.log_multiplier),
                   "calibration_count":int(fitted.calibration_count),"calibration_log_rmse":float(fitted.calibration_log_rmse)},
        "cases":cases,"metrics":metrics,"leave_one_calibration_out":loo,
        "primary_screen":{"threshold":limit,"passed":not failed,"failed_rows":failed},
        "conditional_process":process,
        "timing_s":{"runtime_access_or_initialization":initialization_s,"analysis_including_first_compilation":perf_counter()-started}}


def export(snapshot, result, output):
    """Write a new immutable run directory. Never overwrite a completed run."""
    from calibration.plot import plot_results
    output=Path(output)
    output.mkdir(parents=True,exist_ok=False)
    write_json(output/"inputs.json",snapshot)
    with (output/"predictions.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(result["cases"][0]))
        writer.writeheader()
        writer.writerows(result["cases"])
    plot_results(result,output)
    summary={"schema_version":1,"run_id":output.name,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "protocol_id":snapshot["protocol"]["protocol_id"],"execution_status":"completed",
        "scope":"Experimental support reduced-model calibration and retrospective same-study holdout check. No OpenLB geometry or complete-membrane validation.",
        "units":{"permeability":"m2","particle_diameter":"m","porosity":"1","relative_error":"1",
                 "log_rmse":"1 (natural log)","pressure":"Pa","resistance":"m^-1","shaft_energy":"kWh/m3","timing":"s"},
        "runtime_versions":runtime_versions(),**result,
        "artifact_sha256":{f.name:sha256(f) for f in sorted(output.iterdir())}}
    write_json(output/"summary.json",summary)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,help="New directory; must not already exist")
    parser.add_argument("--require-sources",action="store_true",help="Require every audited local source file; hashes of present sources are always checked")
    args=parser.parse_args()
    snapshot=load_snapshot(require_sources=args.require_sources)
    result=analyze(snapshot)
    output=args.output or ROOT/"results/calibration/runs"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8])
    export(snapshot,result,output)
    print(json.dumps({"output":str(output.resolve()),"execution":"completed", "fit":result["fit"],
        "primary_screen":result["primary_screen"],"holdout":result["metrics"]["primary_holdout"]},indent=2,allow_nan=False))
    # Completion does not imply scientific acceptance. The checker offers
    # --require-screen-pass for workflows that require that separate outcome.


if __name__=="__main__":
    main()
