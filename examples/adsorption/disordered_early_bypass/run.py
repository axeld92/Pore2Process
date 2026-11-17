"""Audit held-out first-contact escape and fixed-conductance clean-film limits."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import crossing
from examples.adsorption.carbon_access_transfer import kernel, run as transfer
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .model import (archive_lags, cdf_at, compose, during_holdout, maximum_cdf_error,
                    plateau_bound, propagation, required_ntu, wall_stages)

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_early_bypass"
COLUMN = transfer.OUTPUT/"runs/20261117T092847Z-6bf8a9d8/summary.json"
COLUMN_SHA = "8a22cc4817aed778218b7d794c5150c5ae0e6b55c5fc27bbce1acf0a9378c32a"
PROBES_SHA = "9a3fe478bcd839c0ba0210a998dcd62bb335978a2b4c71dd886b8479b4bb5167"
FLOWS = ("q233","q625")
STEPS = {"coarse":.002,"fine":.001}
MODES = ("inert","receiver","pair","joint")
LIMIT = .02
TARGET = .05


def inputs():
    if sha(COLUMN)!=COLUMN_SHA:
        raise ValueError("Prior adsorption transfer archive changed")
    record = validate_record(COLUMN)
    base = COLUMN.parent
    plan,frozen = (read_json(base/name) for name in ("plan.json","frozen.json"))
    probes = ROOT/plan["probe_summary"]
    if sha(probes)!=PROBES_SHA or sha(base/"frozen.json")!=record["frozen_sha256"]:
        raise ValueError("Frozen N160 probe or column manifest changed")
    probe = validate_record(probes)
    if any(not probe["cases"][flow]["variance_passed"] for flow in FLOWS):
        raise ValueError("A saved carrier flow lacks the tracer-variance qualification")
    return plan,frozen,probes.parent


def probe_field(probes,flow,part):
    n = read_json(probes/"plan.json")["cases"][flow]["training_count"]
    if part=="training":
        return np.fromfile(probes/flow/"training.f64",dtype="<f8").reshape(n,2)
    if part=="operator":
        return np.fromfile(probes/flow/"operator.f64",dtype="<f8").reshape(2*n,16)
    return np.fromfile(probes/flow/"markers.f64",dtype="<f8",count=16*n).reshape(n,16)


def model_scale(name,frozen):
    file = frozen["cases"][name]["f50_slow01"]["selected"]
    return read_json(COLUMN.parent/file)["result"]["delay_scale"]


def save_pulse(out,key,probabilities):
    for i,p in enumerate(probabilities,1):
        out[f"{key}_{i}"] = p


def run_wall(base,probes):
    jl = runtime()
    jl.include(str(HERE/"states.jl"))
    for flow in FLOWS:
        operator = probe_field(probes,flow,"operator")
        arrivals,hits = np.split(operator,2,axis=1)
        states = np.array(jl.BypassStates.classify(probe_field(probes,flow,"training"),
                          arrivals.copy(order="C")),dtype=int)-1
        stages = wall_stages(arrivals,hits,states)
        data = {}
        for grid,step in STEPS.items():
            for mode in MODES:
                save_pulse(data,f"{grid}_{mode}",propagation(stages[mode],step))
        np.savez_compressed(base/f"{flow}_wall.npz",**data)
        print(flow,"reduced wall-free pulses frozen",flush=True)


def profile(data,grid,mode):
    return [data[f"{grid}_{mode}_{i}"] for i in range(1,9)]


def pulse_summary(pulse,step,holding):
    exit_fraction = float(pulse.sum())
    times = step*np.arange(len(pulse))
    return dict(escape_fraction=exit_fraction,escape_by_three_residence_times=float(cdf_at(pulse,step,3*holding)),
                t05_s=crossing(times,np.cumsum(pulse),TARGET))


def score_wall(base,probes,plan,frozen):
    result = {}
    for flow in FLOWS:
        raw = probe_field(probes,flow,"holdout")
        archive_lags(raw[:,:8])
        condition = plan["cases"][flow+"_l25"]
        p = condition["parameters"]
        blocks,tail = p["blocks"],p["terminal_sections"]
        scale = model_scale(flow+"_l25",frozen)
        holding = p["liquid_volume_m3"]/p["flow_m3_s"]
        ref = during_holdout(raw,STEPS["fine"])
        folded_ref = compose(ref,blocks,tail)
        score = {}
        with np.load(base/f"{flow}_wall.npz") as data:
            for mode in MODES[1:]:
                fine,coarse = [profile(data,k,mode) for k in ("fine","coarse")]
                full = [maximum_cdf_error(a,b,STEPS["fine"],STEPS["fine"])
                        for a,b in zip(fine[2:],ref[2:],strict=True)]
                cuts = np.quantile(raw[:,2:8],[.10,.25],axis=0)
                early = [float(max(abs(cdf_at(a,STEPS["fine"],t)-cdf_at(b,STEPS["fine"],t))))
                         for a,b,t in zip(fine[2:],ref[2:],cuts.T,strict=True)]
                step_errors = [maximum_cdf_error(a,b,STEPS["fine"],STEPS["coarse"])
                               for a,b in zip(fine[2:],coarse[2:],strict=True)]
                composed = compose(fine,blocks,tail)
                score[mode] = dict(escape_by_section=[float(v.sum()) for v in fine],
                    maximum_cdf_error=max(full),maximum_early_error=max(early),cdf_errors=full,early_errors=early,
                    maximum_time_bin_change=max(step_errors),passed=max(full+early)<=LIMIT,
                    column=pulse_summary(composed,STEPS["fine"]*scale,holding),
                    column_cdf_error=maximum_cdf_error(composed,folded_ref,STEPS["fine"],STEPS["fine"]))
            inert = profile(data,"fine","inert")
            np.testing.assert_allclose([v.sum() for v in inert],1.,atol=1e-11,rtol=0.)
        result[flow] = dict(heldout_wall_free_by_section=[float(v.sum()) for v in ref],
            heldout_fastest_tenth_wall_free_dose=[float(cdf_at(v,STEPS["fine"],t)) for v,t in
                zip(ref,np.quantile(raw[:,:8],.1,axis=0),strict=True)],
            physical_time_scale=scale,liquid_residence_s=holding,blocks=blocks,terminal_sections=tail,
            observed_column_projection=pulse_summary(folded_ref,STEPS["fine"]*scale,holding),models=score)
    return result


def film_stages(model,total_ntu,sections,allocate):
    flow = np.array(model.flow)
    contact = np.array(model.area)
    stages = []
    for k,stage in enumerate(model.sections):
        i,j = [np.array(getattr(stage,key),dtype=int)-1 for key in ("donor","receiver")]
        seconds,share = [np.array(getattr(stage,key)) for key in ("seconds","share")]
        # Shares in the saved rules are joint carrier fluxes.
        donor_flux = np.bincount(i,weights=share,minlength=len(flow))
        area = contact[:,k] if allocate=="wall" else flow
        attenuation = 1/(1+total_ntu/sections*area/flow)
        stages.append((i,j,seconds,share/donor_flux[i]*attenuation[j]))
    return stages


def score_film(base,probes,plan,wall):
    result = {}
    dose = {}
    for flow in FLOWS:
        if wall[flow]["observed_column_projection"]["escape_fraction"]>=TARGET:
            continue
        condition = plan["cases"][flow+"_l25"]
        p = condition["parameters"]
        s = condition["scenarios"]["f50_slow01"]["inputs"]
        model = kernel.transport(condition["flow_mL_min"],probes)
        area = 6*p["dry_mass_kg"]/(p["particle_density_kg_m3"]*p["particle_diameter_m"])
        direct = s["direct_fraction"]
        G = area*p["film_m_s"]*(direct+(1-direct)*s["slow_ratio"])
        ntu = G/p["flow_m3_s"]
        sections = condition["sections"]
        bound = plateau_bound(ntu,sections)
        target_ntu = required_ntu(TARGET,sections)
        step = STEPS["fine"]*wall[flow]["physical_time_scale"]
        scores = {}
        for name in ("wall","uniform_per_flow"):
            sample = propagation(film_stages(model,ntu,sections,name),STEPS["fine"])
            pulse = compose(sample,p["blocks"],p["terminal_sections"])
            scores[name] = pulse_summary(pulse,step,wall[flow]["liquid_residence_s"])
            dose[f"{flow}_{name}"] = pulse
        np.testing.assert_allclose(scores["uniform_per_flow"]["escape_fraction"],bound,atol=1e-11,rtol=0.)
        if scores["wall"]["escape_fraction"]<bound-1e-11:
            raise ValueError("Wall-weighted contact violated the fixed-conductance bound")
        result[flow] = dict(total_film_conductance_m3_s=G,total_ntu=ntu,
            best_escape_fraction=bound,minimum_ntu_for_target=target_ntu,
            conductance_increase_required_fraction=max(0.,target_ntu/ntu-1),**scores)
    np.savez_compressed(base/"clean_film.npz",**dose)
    return result


def run():
    plan,frozen,probes = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/name for name in ("run.py","model.py","states.jl","protocol.md")])
    write_json(base/"plan.json",dict(column_sha256=COLUMN_SHA,probes_sha256=PROBES_SHA,steps=STEPS,limit=LIMIT,target=TARGET))
    print(base,flush=True)
    run_wall(base,probes)
    freeze = base/"frozen.json"
    write_json(freeze,{flow:sha(base/f"{flow}_wall.npz") for flow in FLOWS})
    print("Predictions frozen before scoring held-out wall and timing data",sha(freeze),flush=True)
    wall = score_wall(base,probes,plan,frozen)
    film = score_film(base,probes,plan,wall)
    qualified = all(wall[key]["models"]["joint"]["passed"] and key in film
                    and film[key]["best_escape_fraction"]<TARGET for key in FLOWS)
    record = base/"summary.json"
    write_json(record,dict(schema="pore2process.early-bypass.v1",run_id=stamp,
        source_sha256=sources,artifacts_sha256={p.name:sha(p) for p in sorted(base.iterdir())},
        column_sha256=COLUMN_SHA,probes_sha256=PROBES_SHA,frozen_sha256=sha(freeze),
        perfect_sink=wall,clean_film=film,adsorption_rerun_gate_passed=qualified,
        new_adsorption_runs=0,figure7_fits=0,versions=runtime_versions()))
    write_json(OUTPUT/"latest.json",dict(summary=record.relative_to(ROOT).as_posix()))
    print(record,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    plan,frozen,probes = inputs()
    freeze = base/"frozen.json"
    if (sha(freeze)!=record["frozen_sha256"]
            or any(sha(base/f"{flow}_wall.npz")!=digest for flow,digest in read_json(freeze).items())
            or record["perfect_sink"]!=score_wall(base,probes,plan,frozen)):
        raise ValueError("Frozen wall comparison differs on replay")
    with np.load(base/"clean_film.npz") as dose:
        for flow,r in record["clean_film"].items():
            wall = record["perfect_sink"][flow]
            step = STEPS["fine"]*wall["physical_time_scale"]
            for key in ("wall","uniform_per_flow"):
                if r[key]!=pulse_summary(dose[f"{flow}_{key}"],step,wall["liquid_residence_s"]):
                    raise ValueError("Clean-film diagnostics differ on replay")
            sections = plan["cases"][flow+"_l25"]["sections"]
            if r["best_escape_fraction"]!=plateau_bound(r["total_ntu"],sections):
                raise ValueError("Clean-film plateau bound differs on replay")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",type=Path)
    args = p.parse_args()
    check(args.check.resolve()) if args.check else run()
