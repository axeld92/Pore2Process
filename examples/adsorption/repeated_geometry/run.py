"""Run a frozen eight-cell wall-contact prediction and score independent tracers."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.exchange_column.closure import EXCHANGE_SHA, MARKERS
from examples.adsorption.survivor_exposure.run import ARCHIVE, ARCHIVE_SHA
from examples.adsorption.survivor_exposure.column import HOLDOUT
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .predict import LIMITS, PERIODS, comparison, observe, profiles

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/repeated_geometry"
HOLDOUT_SHA = "d98fab904f9584d6953d799f37684a0539642928f75ac03fd6c8cb5121e7e27e"
COUNT = 16384
STOP_S = 12.
PARITY_COUNT = 512


def inputs(base):
    for path, digest in ((HOLDOUT, HOLDOUT_SHA), (ARCHIVE/"summary.json", ARCHIVE_SHA),
                         (MARKERS/"summary.json", EXCHANGE_SHA)):
        if sha(path) != digest:
            raise ValueError("Pinned carrier, wall or two-cell evidence changed")
        validate_record(path)
    holdout = read_json(HOLDOUT)
    if not holdout["scorecard"]["wall_holdout_passed"]:
        raise ValueError("Original wall holdout failed")
    frozen = read_json(HOLDOUT.parent/"frozen_prediction.json")
    wall = read_json(ARCHIVE/"model.json")
    if frozen["dry_carbon_mass_kg"] != wall["dry_carbon_mass_kg"]:
        raise ValueError("Shared dry carbon inventory changed")
    requests = read_json(ARCHIVE/"plan.json")["cases"]
    mesh = read_json(MARKERS/"summary.json")["mesh"]
    period_voxels = round(mesh["length_m"]/2/(mesh["dx_m"]*3))
    period_length = mesh["length_m"]/2
    flow = np.array(frozen["carrier_flow_fractions"])
    predictions = {}
    commands = {}
    for key, case in frozen["cases"].items():
        size = key.split("_")[0]
        fixed = np.array(case["fixed_receiving_mix"])
        np.testing.assert_array_equal(fixed, wall["sizes"][size]["mix"])
        sub = case["models"]["twin_quartiles"]
        groups = sub["subgroups"]
        joint = np.array(sub["carrier_flux"])
        block = joint.reshape(len(flow), groups, len(flow), groups)
        np.testing.assert_allclose(block.sum(axis=(1, 3)), fixed.T*flow[None, :], rtol=2e-14, atol=0.)
        args = (sub["fresh_survival"], sub["flow"], joint, groups)
        repeated = profiles(*args)
        reset = profiles(*args, reset_every=2)
        for name in repeated:
            np.testing.assert_array_equal(repeated[name][:2], reset[name][:2])
        predictions[key] = dict(repeated=repeated, reset_after_two_cells=reset,
                                same_matrix_verified=True, diffusivity_m2_s=case["diffusivity_m2_s"],
                                dt_s=case["dt_s"])
        runtime = requests[key]["command"][:4]
        _, grid, _, dx, d, dt, old_stop, old_count, seed = requests[key]["command"][4:]
        if Path(grid) != MARKERS/"compact" or float(dx) != mesh["dx_m"]*3:
            raise ValueError("Wall run did not use the pinned compact carrier")
        new_seed = int(seed)+3*int(old_count)
        commands[key] = [*runtime, str(HERE/"track.jl"), grid, str(base/(key+".f64")), dx,d,dt,
                         str(STOP_S),str(COUNT),str(new_seed),str(PERIODS),str(period_voxels)]
        if key == "d800_primary":
            parity = [*runtime, str(HERE/"track.jl"), grid, str(base/"parity.f64"), dx,d,dt,
                      old_stop,str(PARITY_COUNT),seed,"2",str(period_voxels)]
    return dict(predictions=predictions, commands=commands, parity_command=parity,
                dry_carbon_mass_kg=frozen["dry_carbon_mass_kg"], period_length_m=period_length,
                carrier_flow_m3_s=mesh["flow_m3_s"], holdout_sha256=HOLDOUT_SHA,
                carrier_sha256=EXCHANGE_SHA, wall_sha256=ARCHIVE_SHA, periods=PERIODS,
                marker_count=COUNT, limits=LIMITS)


def parity_check(base):
    new = np.fromfile(base/"parity.f64", dtype="<f8").reshape(PARITY_COUNT, 4)
    old = np.fromfile(ARCHIVE/"d800_primary.f64", dtype="<f8").reshape(-1, 4)[:PARITY_COUNT]
    new[:, (1, 3)] -= new[:, [0]]
    np.testing.assert_array_equal(new, old)


def score(base, plan):
    parity_check(base)
    result = {}
    for key, frozen in plan["predictions"].items():
        data = np.fromfile(base/(key+".f64"), dtype="<f8").reshape(COUNT, 2*PERIODS)
        tracked = observe(data)
        recovered = tracked["recovery"][-1] >= .999 and tracked["maximum_unresolved_fraction"] <= .001
        result[key] = dict(observed=tracked, inert_recovery_passed=recovered,
                           repeated=comparison(tracked, frozen["repeated"]),
                           reset=comparison(tracked, frozen["reset_after_two_cells"]))
    refinement = {size: comparison(result[size+"_primary"]["observed"], result[size+"_half_step"]["observed"])
                  for size in ("d800", "d450")}
    passed = all(case["inert_recovery_passed"] and case["repeated"]["passed"] for case in result.values())
    return dict(cases=result, time_step=refinement, longer_holdout_passed=passed and
                all(case["passed"] for case in refinement.values()), archived_two_cell_parity_markers=PARITY_COUNT)


def run():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    plan = inputs(base)
    sources = source_manifest([HERE/name for name in ("run.py","track.jl","predict.py","protocol.md","test_track.jl")]
                              +[HERE.parent/name for name in ("carbon_contact/track_contacts.jl","path_exchange/track.jl")])
    for name in sources:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name,target)
    frozen = base/"frozen_prediction.json"
    write_json(frozen, plan)
    digest = sha(frozen)
    print("Long tracer prediction frozen", digest, base, flush=True)
    for key, cmd in {"parity":plan["parity_command"], **plan["commands"]}.items():
        timing = execute(cmd,base/(key+".log"),timeout=1800)
        write_json(base/(key+"_execution.json"),timing)
        if key == "parity":
            parity_check(base)
        print(key,timing,flush=True)
    result = score(base, plan)
    path = base/"summary.json"
    write_json(path, dict(schema="pore2process.repeated-geometry.v1",run_id=stamp,
                          source_sha256=sources,frozen_prediction_sha256=digest,
                          artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
                          analysis=result,new_openlb_runs=0,new_column_runs=0,figure7_fits=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Complete",path,"holdout passed",result["longer_holdout_passed"],flush=True)


def check(path):
    saved = validate_record(path)
    base = path.parent
    frozen = base/"frozen_prediction.json"
    plan = inputs(base)
    if sha(frozen) != saved["frozen_prediction_sha256"] or read_json(frozen) != plan:
        raise ValueError("Frozen prediction differs on replay")
    if score(base,plan) != saved["analysis"]:
        raise ValueError("Holdout score differs on replay")
    print(dict(sources=len(saved["source_sha256"]),artifacts=len(saved["artifacts_sha256"]),
               passed=saved["analysis"]["longer_holdout_passed"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
