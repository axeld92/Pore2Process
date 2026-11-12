"""Freeze first-period wall-exposure models and check the saved second period."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.contact_column.inputs import MARKERS, rank_bins
from examples.adsorption.transport_accuracy.run import validate_record
from examples.adsorption.wall_uptake.estimate import SHARED, SHARED_SHA
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .predict import LIMITS, QUARTILES, predict, score, states

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/survivor_exposure"
ARCHIVE = ROOT/"results/adsorption/wall_uptake/runs/20261112T143948Z-126af769"
ARCHIVE_SHA = "d2cfb731915c913660317a1c58b51ef4256067ce2bf230f59491082fe5b47d79"
GROUPS = {"parent_only": 1, "twin_quartiles": QUARTILES}


def inputs():
    for path, expected in ((ARCHIVE, ARCHIVE_SHA), (SHARED, SHARED_SHA)):
        if sha(path/"summary.json") != expected:
            raise ValueError("Pinned wall-uptake or shared-carbon cohort changed")
        validate_record(path/"summary.json")
    marker = validate_record(MARKERS/"summary.json")
    twins = np.fromfile(MARKERS/"advection.f64", dtype="<f8").reshape(-1, 2)
    if not np.isfinite(twins).all() or np.any(twins <= 0.):
        raise ValueError("Invalid advection-only twin times")
    plan = read_json(ARCHIVE/"plan.json")
    model = read_json(ARCHIVE/"model.json")
    shared = read_json(SHARED/"model.json")
    if model["dry_carbon_mass_kg"] != shared["dry_carbon_mass_kg"]:
        raise ValueError("Dry carbon inventory differs from the shared reference")
    flow = np.array(model["closure"]["flow_fractions"])
    cases = {}
    for key, spec in plan["cases"].items():
        size = key.split("_")[0]
        setup = model["sizes"][size]
        if (setup["mix"] != shared["sizes"][size]["mix"]
                or setup["parameters"] != shared["sizes"][size]["parameters"]
                or spec["command"][-2:] != marker["requests"]["advection"][-2:]):
            raise ValueError("Péclet exchange, particle inputs or marker identity changed")
        cases[key] = (spec, np.array(setup["mix"]))
    return flow, twins, model["dry_carbon_mass_kg"], cases


def fields(key, count):
    return np.memmap(ARCHIVE/(key+".f64"), dtype="<f8", mode="r", shape=(count, 4))


def freeze(flow, twins, mass, cases):
    result = {}
    for key, (spec, receiving) in cases.items():
        data = fields(key, len(twins))
        transit = np.array(data[:, :2])
        hit = np.array(data[:, 2])
        first_survives = np.isnan(hit)
        if (not np.isfinite(transit).all() or np.any(transit <= 0.) or np.isinf(hit).any()
                or np.any(hit[~first_survives] < 0.)
                or np.any(hit[~first_survives] > transit[~first_survives, 0]+1e-12)):
            raise ValueError("Invalid transit or first-cell hit time")
        models = {name: predict(states(transit, twins, flow, groups),
                                first_survives, flow, receiving, groups)
                  for name, groups in GROUPS.items()}
        result[key] = dict(dt_s=spec["dt_s"], diffusivity_m2_s=spec["diffusivity_m2_s"],
                           data_sha256=sha(ARCHIVE/(key+".f64")),
                           fixed_receiving_mix=receiving.tolist(), models=models)
    return dict(dry_carbon_mass_kg=mass, carrier_flow_fractions=flow.tolist(),
                marker_count=len(twins), parent_wall_uptake_sha256=ARCHIVE_SHA,
                parent_shared_carbon_sha256=SHARED_SHA, cases=result)


def evaluate(frozen, flow, twins):
    result = {}
    advection = rank_bins(twins, np.r_[0., np.cumsum(flow)])[:, 1]
    old = read_json(ARCHIVE/"model.json")
    for key, case in frozen["cases"].items():
        if sha(ARCHIVE/(key+".f64")) != case["data_sha256"]:
            raise ValueError("Marker archive changed after prediction freeze")
        data = fields(key, len(twins))
        first_survives = np.isnan(data[:, 2])
        hit = data[:, 3]
        second_hit = np.isfinite(hit)
        if (np.isinf(hit).any() or np.any(hit[second_hit] < 0.)
                or np.any(hit[second_hit] > data[second_hit, 1]+1e-12)):
            raise ValueError("Invalid held-out hit time")
        size, variant = key.split("_", 1)
        prior = old["sizes"][size]["uptake"][variant]
        scores = {}
        for name, model in case["models"].items():
            labels = states(data[:, :2], twins, flow, model["subgroups"])
            scores[name] = score(model, labels, advection, first_survives, second_hit, flow)
        result[key] = dict(models=scores, original_constant_film_errors=prior["errors"])
    passed = all(case["models"]["twin_quartiles"]["passed"] for case in result.values())
    return dict(cases=result, limits=LIMITS, wall_holdout_passed=passed)


def run():
    flow, twins, mass, cases = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    sources = source_manifest([HERE/name for name in ("predict.py", "protocol.md", "run.py")])
    for name in sources:
        copy = base/"source_snapshot"/name
        copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, copy)
    frozen_path = base/"frozen_prediction.json"
    frozen = freeze(flow, twins, mass, cases)
    write_json(frozen_path, frozen)
    frozen_sha = sha(frozen_path)
    print("Prediction frozen before scoring the second-cell hits:", frozen_sha, flush=True)
    scorecard = evaluate(frozen, flow, twins)
    summary = base/"summary.json"
    write_json(summary, dict(schema="pore2process.survivor-exposure.v1", run_id=stamp,
                             source_sha256=sources, frozen_prediction_sha256=frozen_sha,
                             artifacts_sha256={str(p.relative_to(base)): sha(p)
                                               for p in base.rglob("*") if p.is_file()},
                             scorecard=scorecard, figure7_reruns=0))
    write_json(OUTPUT/"latest.json", dict(summary=str(summary.relative_to(ROOT))))
    for key, case in scorecard["cases"].items():
        print(key, case["models"]["twin_quartiles"]["errors"], flush=True)
    print("Wall holdout passed:", scorecard["wall_holdout_passed"], summary, flush=True)


def check(path):
    record = validate_record(path)
    flow, twins, mass, cases = inputs()
    frozen_path = path.parent/"frozen_prediction.json"
    frozen = read_json(frozen_path)
    if sha(frozen_path) != record["frozen_prediction_sha256"] or frozen != freeze(flow, twins, mass, cases):
        raise ValueError("Frozen prediction differs on replay")
    if record["scorecard"] != evaluate(frozen, flow, twins):
        raise ValueError("Holdout score differs on replay")
    print(dict(sources=len(record["source_sha256"]), artifacts=len(record["artifacts_sha256"]),
               wall_holdout_passed=record["scorecard"]["wall_holdout_passed"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
