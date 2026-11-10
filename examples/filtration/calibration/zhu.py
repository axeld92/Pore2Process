"""Real Julia fit/holdout test of the digitized Zhu Figure 6 M1 series."""
from pore2process_demo.paths import project_path
from pore2process_demo.provenance import implementation_files, source_manifest
import argparse
import csv
from datetime import datetime, timezone
from itertools import product
import math
from pathlib import Path
from time import perf_counter
from uuid import uuid4

from examples.filtration.calibration.run import ROOT, read_json, write_json, sha256
from examples.filtration.calibration.check import close_tree
from examples.filtration.calibration.digitize_zhu import pixels_to_flux
from pore2process_demo.bridge import runtime, runtime_versions


def load_snapshot(require_sources=False,source_reading="jpeg"):
    if source_reading not in ("jpeg","pdf"):
        raise ValueError("Source reading must be jpeg or pdf.")
    source=read_json(ROOT/'literature/sources/zhu2021/extraction.json')
    status=[]
    for name,digest in source["source"]["files_sha256"].items():
        p=project_path(ROOT/source["source"]["folder"]/name)
        if p.exists():
            if sha256(p)!=digest:
                raise ValueError(f"Source changed; review first: {p}")
            state="sha256 verified"
        elif require_sources:
            raise ValueError(f"Required source missing: {p}")
        else:
            state="missing; retained audited extraction used"
        status.append({"path":str(p.relative_to(ROOT)),"sha256":digest,"status":state})
    protocol_path="examples/filtration/calibration/zhu_protocol_pdf.json" if source_reading=="pdf" else "examples/filtration/calibration/zhu_protocol.json"
    data_path='literature/sources/zhu2021/figure6_pdf_points.json' if source_reading=="pdf" else 'literature/sources/zhu2021/figure6_jpeg_points.json'
    names=["examples/filtration/calibration/zhu.py",protocol_path,"examples/filtration/calibration/digitize_zhu.py",
           'literature/sources/zhu2021/extraction.json',data_path,
           "julia/src/SeriesCalibration.jl","julia/src/Pore2ProcessDemo.jl",
           "src/pore2process_demo/bridge.py","uv.lock","julia/Project.toml","julia/Manifest.toml"]
    names = implementation_files(names)
    return {"protocol":read_json(project_path(ROOT/protocol_path)),"source":source,"reading_source":source_reading,
            "digitized":read_json(project_path(ROOT/data_path)),"source_files":status,
            "implementation_sha256":source_manifest(names)}


def validate(snapshot):
    data=snapshot["digitized"]
    points=data["points"]
    config_key="digitization_pdf" if snapshot.get("reading_source","jpeg")=="pdf" else "digitization"
    if data["config"]!=snapshot["source"][config_key]:
        raise ValueError("Digitization configuration differs from audited source.")
    if len(points)!=data["config"]["expected_point_count"]:
        raise ValueError("Digitized point count differs from audit.")
    last=-math.inf
    for i,row in enumerate(points,1):
        if row["point_id"]!=i or row["pixel_x"]<=last:
            raise ValueError("Point IDs must follow strictly increasing support-flux rank.")
        last=row["pixel_x"]
        js,jm=pixels_to_flux(data["config"],row["pixel_x"],row["pixel_y"])
        if not math.isclose(js,row["support_flux_m_h"],rel_tol=1e-13) or not math.isclose(jm,row["total_flux_LMH"],rel_tol=1e-13):
            raise ValueError("Digitized flux does not match saved pixel/axis coordinates.")
        if not 0<jm<1000*js:
            raise ValueError("Paired total flux must be positive and smaller than support flux.")


def metrics(rows):
    errors=[r["relative_flux_error"] for r in rows]
    return {"count":len(rows),"mean_absolute_relative_error":sum(map(abs,errors))/len(rows),
            "max_absolute_relative_error":max(map(abs,errors)),"mean_relative_error":sum(errors)/len(rows)}


def analyze(snapshot,with_sensitivity=True):
    validate(snapshot)
    p=snapshot["protocol"]
    config=snapshot["digitized"]["config"]
    points=snapshot["digitized"]["points"]
    dp=snapshot["source"]["measurement_conditions"]["test_pressure_drop_Pa"]
    started=perf_counter()
    jl=runtime()
    initialization_s=perf_counter()-started
    started=perf_counter()
    vector=jl.seval("Vector{Float64}")
    train=[r for r in points if r["point_id"]%2]
    def fit_fluxes(fluxes):
        return jl.series_impedance_fit(vector([js/3600/dp for js,jm in fluxes]),
                                       vector([jm/3.6e6/dp for js,jm in fluxes]))
    def fit_points(selected):
        return fit_fluxes([(r["support_flux_m_h"],r["total_flux_LMH"]) for r in selected])
    fitted=fit_points(train)
    Z=float(fitted.fixed_impedance_Pa_s_m)
    limiting=dp*float(fitted.limiting_permeance_m_Pa_s)
    cases=[]
    for row in points:
        ps=row["support_flux_m_h"]/3600/dp
        measured=row["total_flux_LMH"]/3.6e6
        value=jl.series_flux(ps,Z,dp)
        prediction=float(value.flux_m_s)
        cases.append({**row,"role":'calibration' if row["point_id"]%2 else "holdout",
            "support_permeance_m_Pa_s":ps,"measured_flux_m_s":measured,
            "predicted_flux_m_s":prediction,"predicted_flux_LMH":prediction*3.6e6,
            "relative_flux_error":prediction/measured-1,
            "non_support_fraction":float(value.non_support_fraction)})
    sensitivity=None
    if with_sensitivity:
        # Bounds on each pixel coordinate and coherent bounds on four axis
        # anchors. Monotonicity gives min/max training intercepts at corners.
        delta=p["digitization"]["marker_center_bound_px"]
        anchor_delta=p["digitization"]["axis_anchor_bound_px"]
        Zs=[]
        largest=0.
        for shifts in product((-anchor_delta,anchor_delta),repeat=4):
            anchors={key:val+shift for (key,val),shift in zip(config["axis_anchors"].items(),shifts,strict=True)}
            extremes=[]
            for sign in (-1,1):
                fluxes=[pixels_to_flux(config,r["pixel_x"]+sign*delta,r["pixel_y"]+sign*delta,anchors) for r in train]
                value=float(fit_fluxes(fluxes).fixed_impedance_Pa_s_m)
                Zs.append(value); extremes.append(value)
            for row in cases:
                if row["role"]!="holdout":
                    continue
                for dx,dy in product((-delta,delta),repeat=2):
                    js,jm=pixels_to_flux(config,row["pixel_x"]+dx,row["pixel_y"]+dy,anchors)
                    for z in extremes:
                        pred=float(jl.series_flux(js/3600/dp,z,dp).flux_m_s)*3.6e6
                        largest=max(largest,abs(pred/jm-1))
        sensitivity={"fixed_impedance_range_Pa_s_m":[min(Zs),max(Zs)],
                     "worst_corner_holdout_absolute_relative_error":largest,
                     "scope":"Deterministic pixel/axis reading sensitivity; not measurement uncertainty or statistical confidence."}
    loads=[]
    settings=p["conditional_process"]
    js_min,js_max=settings["support_flux_range_m_h"]
    support_targets=[js_min,9*limiting*3600,30.,js_max]
    for js in sorted(set(support_targets)):
        v=jl.impedance_load(js/3600/dp,Z,settings["filtration_flux_m_s"],settings["membrane_area_m2"],settings["pump_efficiency_fraction"])
        fields=("pressure_drop_Pa","flow_m3_s","shaft_power_W","specific_shaft_energy_kWh_m3","non_support_fraction","support_fraction")
        loads.append({"support_flux_at_1bar_m_h":js,**{name:float(getattr(v,name)) for name in fields},
                      "above_authors_30_m_h_strength_screen":js>30})
    baseline=loads[0]["specific_shaft_energy_kWh_m3"]
    for load in loads:
        load["energy_saving_vs_low_support_fraction"]=1-load["specific_shaft_energy_kWh_m3"]/baseline
    groups={role:metrics([r for r in cases if r["role"]==role]) for role in ('calibration',"holdout")}
    limit=p["primary_screen"]["max_absolute_relative_flux_error"]
    failures=[c["point_id"] for c in cases if c["role"]=="holdout" and abs(c["relative_flux_error"])>limit]
    return {"fit":{"fixed_impedance_Pa_s_m":Z,"calibration_count":len(train),"limiting_flux_LMH":limiting*3.6e6,
                   "impedance_rms_residual_Pa_s_m":float(fitted.impedance_rms_residual_Pa_s_m)},
        "cases":cases,"metrics":groups,"digitization_sensitivity":sensitivity,
        "screen":{"nominal_passed":not failures,"failed_holdout_ids":failures,"relative_error_limit":limit,
                  "all_reading_corners_passed":None if sensitivity is None else sensitivity["worst_corner_holdout_absolute_relative_error"]<=limit},
        "threshold_support_flux_at_1bar_m_h":{"non_support_50_percent":limiting*3600,"non_support_90_percent":9*limiting*3600},
        "conditional_process":loads,
        "reported_endpoint_energy_translation":{
            name:{"reported_fixed_pressure_flux_increase_fraction":gain,
                  "conditional_fixed_flux_energy_saving_fraction":1-1/(1+gain)}
            for name,gain in snapshot["source"]["reported_endpoint_flux_increase_fractions"].items()},
        "timing_s":{"runtime_access_or_initialization":initialization_s,"analysis_including_first_compilation":perf_counter()-started}}


def plot(result,output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout="constrained")
    ax,b=axes
    for role,marker,color in (('calibration',"o","#16769a"),("holdout","s","#cf6712")):
        rows=[c for c in result["cases"] if c["role"]==role]
        ax.scatter([r["support_flux_m_h"] for r in rows],[r["total_flux_LMH"] for r in rows],marker=marker,c=color,label=role.capitalize(),s=40)
    rows=sorted(result["cases"],key=lambda r:r["support_flux_m_h"])
    ax.plot([r["support_flux_m_h"] for r in rows],[r["predicted_flux_LMH"] for r in rows],color="black",lw=1,label="Julia series-law prediction")
    ax.axhline(result["fit"]["limiting_flux_LMH"],color="#777777",ls="--",label="Fitted limit")
    ax.set(xlabel="Support flux at 1 bar (m³ m⁻² h⁻¹)",ylabel="Whole-membrane flux at 1 bar (L m⁻² h⁻¹)",title="Zhu Figure 6 M1: digitized data")
    ax.legend(fontsize=8)
    held=[r for r in rows if r["role"]=="holdout"]
    b.axhspan(-10,10,color="#dddddd",alpha=.5,label="±10% project screen")
    b.scatter([r["support_flux_m_h"] for r in held],[100*r["relative_flux_error"] for r in held],marker="s",color="#cf6712")
    b.axhline(0,color="black",lw=1)
    b.set(xlabel="Support flux at 1 bar (m³ m⁻² h⁻¹)",ylabel="(Predicted / digitized flux − 1) (%)",title="Held-out points; no experimental error bars")
    b.legend(fontsize=8)
    fig.savefig(output/"heldout_validation.png",dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.3),layout="constrained")
    loads=result["conditional_process"]
    x=[r["support_flux_at_1bar_m_h"] for r in loads]
    axes[0].plot(x,[100*r["energy_saving_vs_low_support_fraction"] for r in loads],"o-")
    axes[0].set(xscale="log",xlabel="Support flux at 1 bar (m³ m⁻² h⁻¹)",ylabel="Conditional shaft-energy saving (%)",title="Same fitted non-support layer, fixed 100 LMH")
    axes[1].plot(x,[100*r["non_support_fraction"] for r in loads],"o-")
    axes[1].axhline(90,color="#777777",ls="--")
    axes[1].set(xscale="log",xlabel="Support flux at 1 bar (m³ m⁻² h⁻¹)",ylabel="Effective non-support share (%)",title="Share includes possible interface/access effects")
    for a in axes:
        a.grid(alpha=.2)
    fig.savefig(output/"conditional_energy.png",dpi=180)
    plt.close(fig)


def export(snapshot,result,directory):
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=False)
    write_json(directory/"inputs.json",snapshot)
    with (directory/"predictions.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(result["cases"][0]))
        writer.writeheader(); writer.writerows(result["cases"])
    plot(result,directory)
    summary={"schema_version":1,"run_id":directory.name,"timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "execution_status":"completed","protocol_id":snapshot["protocol"]["protocol_id"],
        "units":{"permeance":"m/(Pa s)","impedance":"Pa s/m","flux":"m/s at Julia boundary; m/h or LMH explicitly named in display columns","pressure":"Pa","shaft_power":"W","specific_shaft_energy":"kWh/m3","fractions":"1"},
        "scope":snapshot["protocol"]["scientific_scope"],"runtime_versions":runtime_versions(),
        **result,"artifact_sha256":{f.name:sha256(f) for f in sorted(directory.iterdir())}}
    write_json(directory/"summary.json",summary)
    return summary


def check(directory):
    directory=Path(directory)
    saved=read_json(directory/"summary.json")
    required={"inputs.json","predictions.csv","heldout_validation.png","conditional_energy.png"}
    if set(saved["artifact_sha256"])!=required:
        raise ValueError("Missing or unexpected artifact declarations")
    for name,digest in saved["artifact_sha256"].items():
        if sha256(directory/name)!=digest:
            raise ValueError(f"Artifact hash mismatch: {name}")
    snapshot=read_json(directory/"inputs.json")
    if saved["protocol_id"]!=snapshot["protocol"]["protocol_id"]:
        raise ValueError("Summary protocol differs from saved input")
    replay=analyze(snapshot)
    for key in ("fit","cases","metrics","digitization_sensitivity","screen","threshold_support_flux_at_1bar_m_h","conditional_process","reported_endpoint_energy_translation"):
        close_tree(saved[key],replay[key],key)
    with (directory/"predictions.csv").open(newline="") as f:
        records=list(csv.DictReader(f))
    if len(records)!=len(replay["cases"]):
        raise ValueError("CSV row count differs")
    for row,expected in zip(records,replay["cases"],strict=True):
        for key,value in expected.items():
            actual=None if value is None and row[key]=="" else float(row[key]) if isinstance(value,(float,int)) else row[key]
            close_tree(actual,value,f"CSV.{key}")
    return {"integrity_and_numerical_replay":"passed","screen":replay["screen"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--check",type=Path)
    parser.add_argument("--require-sources",action="store_true")
    parser.add_argument("--source",choices=("jpeg","pdf"),default="pdf",help="PDF native image by default; jpeg reproduces the earlier reading")
    args=parser.parse_args()
    if args.check:
        print(check(args.check)); return
    snapshot=load_snapshot(args.require_sources,args.source)
    result=analyze(snapshot)
    output=args.output or ROOT/"results/zhu_calibration/runs"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8])
    export(snapshot,result,output)
    print({"output":str(output),"fit":result["fit"],"holdout":result["metrics"]["holdout"],"screen":result["screen"],"digitization_sensitivity":result["digitization_sensitivity"]})


if __name__=="__main__":
    main()
