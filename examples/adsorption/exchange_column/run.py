"""Run conservative exchange-column predictions before reading Figure 7 markers."""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
import math
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare, crossing
from examples.adsorption.film_particle.bridge import call, module as particle_module
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .closure import closure

HERE = Path(__file__).parent
OUTPUT = ROOT / "results/adsorption/exchange_column"
JULIA = ROOT / "julia/adsorption"
PARTICLE = ROOT / "examples/adsorption/film_particle/reference/film_particle_reference.json"
SOURCE = ROOT / "literature/sources/alhamed2009/film_particle_extraction.json"
PIXELS = SOURCE.with_name("figure7_pixels.csv")
BATCH_SHA = "7fdc697c9d0b487fe50b914005c914ce8ed06720fe9dab08b7c58770ef1b8bf2"
SOURCE_SHA = "14e362ceda76349c80fbdad1f1a8ebb4298a1adf7d6d4d56e798afae3df3fc5c"
CASES = {"d800": ("calibration", .0008), "d450": ("holdout", .00045)}
LEVELS = {"d800": [(2, 96, 4., .04), (4, 128, 2., .02)],
          "d450": [(1, 24, 4., .08), (2, 48, 2., .04), (4, 96, 1., .02)]}


@lru_cache(maxsize=1)
def solver():
    particle_module()
    jl = runtime()
    jl.include(str(JULIA / "ExchangeColumn.jl"))
    jl.seval("ExchangeColumn.BLAS.set_num_threads(1)")
    return jl.ExchangeColumn.solve_column


def source_inputs():
    if sha(SOURCE) != SOURCE_SHA or sha(PARTICLE) != BATCH_SHA:
        raise ValueError("Frozen source or batch qualification changed")
    source = read_json(SOURCE)
    batch = read_json(PARTICLE)
    checks = batch["physical_status"]
    if not checks["source_equilibrium_audit_passed"] or not checks["batch_particle_size_comparison_passed"]:
        raise ValueError("Frozen particle model is unqualified")
    f = source["figure_7"]
    return source, f, source["source_parameters"], batch["batch_fit"]


def parameters(cl, diameter, level, f, sp, batch):
    L, A, Q = (f[k] for k in ("length_m", "area_m2", "flow_m3_s"))
    eps = sp["porosity"]
    u = Q / A
    length = cl["reference_period_length_m"] * diameter / .0008
    periods = round(L / length)
    film = call("film_coefficient", diameter_m=diameter, superficial_velocity_m_s=u,
                viscosity_Pa_s=sp["viscosity_Pa_s"], density_kg_m3=sp["density_kg_m3"],
                molecular_diffusivity_m2_s=sp["molecular_diffusivity_m2_s"], porosity=eps,
                interstitial=False, multiplier=1.)
    subcells, radial, dt, startup_dt = level
    p = dict(periods=periods, subcells=subcells, length_m=L, area_m2=A, porosity=eps,
             superficial_velocity_m_s=u, inlet_concentration_kg_m3=f["feed_kg_m3"],
             particle_diameter_m=diameter, particle_density_kg_m3=sp["particle_density_kg_m3"],
             diffusivity_m2_s=batch["diffusivity_m2_s"], film_m_s=film["film_m_s"],
             capacity_kg_kg=sp["capacity_kg_kg"], affinity_m3_kg=sp["affinity_m3_kg"],
             radial_cells=radial, t_end_s=18000., dt_s=dt, startup_s=24., startup_dt_s=startup_dt)
    reference_Pe = cl["reference_period_length_m"]**2 / (cl["reference_liquid_transit_s"]*sp["molecular_diffusivity_m2_s"])
    trace = dict(period_length_m=L/periods, scaled_fcc_period_length_m=length,
                 period_nominal_transit_s=eps*L/(u*periods),
                 period_Peclet=u/eps*(L/periods)/sp["molecular_diffusivity_m2_s"], reference_period_Peclet=reference_Pe,
                 column_film=film, liquid_volume_m3=eps*A*L, dry_inventory_kg=f["derived_dry_mass_kg"],
                 bulk_density_kg_m3=sp["bulk_density_kg_m3"], batch_film_ref_m_s=batch["batch_film_ref_m_s"])
    return p, trace


def solve(cl, p, exchange):
    counts = np.asarray(cl["transition_counts"], dtype=np.float64)
    mix = counts.T / counts.sum(axis=0)[:,None] if exchange else np.eye(len(counts))
    jl = runtime()
    result = solver()(np.asarray(cl["flow_fractions"]), np.asarray(cl["normalized_transit_times"]), mix, **p)
    out = {}
    for field in jl.fieldnames(jl.typeof(result)):
        key = str(field)
        value = getattr(result, key)
        out[key] = np.array(value, dtype=np.float64, copy=True).tolist() if bool(jl.isa(value,jl.AbstractArray)) else float(value)
    return out


def validate_solution(result, p, trace):
    low, high = result["concentration_range_kg_m3"]
    q_low, q_high = result["loading_range_kg_kg"]
    if (not all(math.isfinite(x) for x in (low, high, q_low, q_high, result["normalized_mass_residual"],
                                          result["normalized_interface_flux_mismatch"]))
            or low < -1e-12 or high > p["inlet_concentration_kg_m3"] + 1e-10
            or q_low < -1e-12 or q_high > p["capacity_kg_kg"] + 1e-10
            or result["normalized_mass_residual"] > 1e-8
            or result["normalized_interface_flux_mismatch"] > 1e-8):
        raise ValueError("Column conservation or physical bounds failed")
    np.testing.assert_allclose([result["physical_adsorbent_mass_kg"],result["liquid_volume_m3"]],
                               [trace["dry_inventory_kg"],trace["liquid_volume_m3"]],rtol=1e-12,atol=0)


def refine(coarse, fine, feed):
    t = np.asarray(coarse["time_s"])
    c = np.asarray(coarse["outlet_concentration_kg_m3"]) / feed
    tf = np.asarray(fine["time_s"])
    cf = np.asarray(fine["outlet_concentration_kg_m3"]) / feed
    errors = dict(max_fraction=float(max(abs(c - np.interp(t, tf, cf)))))
    for threshold in (.05,.5):
        a, b = crossing(t,c,threshold),crossing(tf,cf,threshold)
        errors[f"t{round(100*threshold):02}"] = abs(a/b - 1) if a is not None and b is not None else None
    passed = errors["max_fraction"] <= .005 and all(errors[k] is not None and errors[k] <= .02 for k in ("t05","t50"))
    return dict(errors=errors, passed=passed)


def score(base, frozen, pixels_hash):
    if sha(PIXELS) != pixels_hash:
        raise ValueError("The Figure 7 digitization changed")
    out = {}
    for size, (series,_) in CASES.items():
        names = frozen["cases"][size]
        a = read_json(base / names["primary"])
        b = read_json(base / names["control"])
        out[size] = dict(exchange=compare(a["result"],series), persistent=compare(b["result"],series),
                         refinement=names["refinement"])
        out[size]["passed"] = names["refinement"][-1]["passed"] and column_passed(out[size]["exchange"])
    return out


def run():
    cl = closure()
    _, f, sp, batch = source_inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    base = OUTPUT / "runs" / stamp
    base.mkdir(parents=True)
    files = [HERE/name for name in ("run.py","closure.py","protocol.md")]
    files += [JULIA/name for name in ("FilmParticle.jl","ExchangeColumn.jl","test_exchange_column.jl")]
    files += [PARTICLE,SOURCE,HERE.parent/"alhamed/source.py",HERE.parent/"alhamed/study.py",
              HERE.parent/"path_exchange/predict.py",HERE.parent/"residence_time_holdout/paths.py"]
    manifest = source_manifest(files)
    for name in manifest:
        target = base / "source_snapshot" / name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT / name,target)
    write_json(base/"closure.json",cl)
    write_json(base/"plan.json",dict(column_parameters_fitted=0,source_sha256=manifest,levels=LEVELS))
    print("Exchange column",base,"classes",len(cl["flow_fractions"]),flush=True)
    frozen = dict(cases={}, hashes={}, closure_sha256=sha(base/"closure.json"))
    def save(name,p,trace,exchange):
        result = solve(cl,p,exchange)
        validate_solution(result,p,trace)
        filename = name + ".json"
        write_json(base/filename,dict(inputs=p,trace=trace,exchange=exchange,result=result))
        frozen["hashes"][filename] = sha(base/filename)
        print(name,"seconds",result["total_s"],"mass",result["normalized_mass_residual"],flush=True)
        return result,filename
    for size,(_,diameter) in CASES.items():
        levels = []
        prior = None
        for index, level in enumerate(LEVELS[size]):
            p, trace = parameters(cl,diameter,level,f,sp,batch)
            result, filename = save(f"{size}_exchange_level{index}",p,trace,True)
            if prior is not None:
                levels.append(refine(prior,result,f["feed_kg_m3"]))
                print(size,"refinement",levels[-1],flush=True)
                if levels[-1]["passed"]:
                    break
            prior = result
        _, control = save(f"{size}_persistent",p,trace,False)
        frozen["cases"][size] = dict(primary=filename,control=control,refinement=levels)
    write_json(base/"frozen_predictions.json",frozen)
    frozen_sha = sha(base/"frozen_predictions.json")
    print("Frozen before Figure 7 scoring",frozen_sha,flush=True)
    comparisons = score(base,frozen,f["pixel_csv_sha256"])
    record = dict(schema="pore2process.exchange-column.v1",run_id=stamp,source_sha256=manifest,
                  artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
                  versions=runtime_versions(),frozen_prediction_sha256=frozen_sha,source_pixel_sha256=f["pixel_csv_sha256"],
                  figure7=comparisons,qualified=all(c["passed"] for c in comparisons.values()),
                  new_column_fitted_parameters=0,new_OpenLB_runs=0,new_DWSIM_runs=0)
    path = base/"summary.json"
    write_json(path,record)
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Completed",path,comparisons,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    frozen = read_json(base/"frozen_predictions.json")
    if sha(base/"frozen_predictions.json") != record["frozen_prediction_sha256"] or read_json(base/"closure.json") != closure():
        raise ValueError("Frozen predictions or source exchange closure changed")
    for name,digest in frozen["hashes"].items():
        if sha(base/name) != digest:
            raise ValueError("A frozen column prediction changed")
        case = read_json(base/name)
        validate_solution(case["result"],case["inputs"],case["trace"])
    if score(base,frozen,record["source_pixel_sha256"]) != record["figure7"]:
        raise ValueError("Figure 7 comparison differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"]),
               qualified=record["qualified"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args=parser.parse_args()
    check(args.check.resolve()) if args.check else run()
