"""Summarize continued fast transit and wall avoidance in the saved holdout."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .predict import FAST, PERIODS
from .run import COUNT, HERE, OUTPUT


def path_persistence(data):
    arrivals, hits = np.split(data, 2, axis=1)
    transit = np.column_stack((arrivals[:, 0], np.diff(arrivals, axis=1)))
    order = np.argsort(transit, axis=0, kind="stable")
    fast = np.zeros(transit.shape, dtype=bool)
    np.put_along_axis(fast, order[:round(FAST*len(data))], True, axis=0)
    # Only completed local transits can contribute to persistence.
    fast &= np.isfinite(transit)
    persistent = fast.all(axis=1)
    wall_free = np.isfinite(arrivals[:, -1]) & ~np.isfinite(hits).any(axis=1)
    slow = ~fast[:, 0] & np.isfinite(arrivals[:, 0])
    return dict(initially_fast_count=int(fast[:, 0].sum()),
                persistent_fast_count=int(persistent.sum()),
                persistent_fraction_of_initially_fast=float(persistent.sum()/fast[:, 0].sum()),
                persistent_wall_free_fraction=float(wall_free[persistent].mean()),
                initially_slow_wall_free_fraction=float(wall_free[slow].mean()))


def summarize(path):
    record = validate_record(path)
    plan = read_json(path.parent/"frozen_prediction.json")
    analysis = record["analysis"]
    cases = {}
    for key, case in analysis["cases"].items():
        data = np.fromfile(path.parent/(key+".f64"), dtype="<f8").reshape(COUNT, 2*PERIODS)
        cases[key] = dict(observed=case["observed"], persistence=path_persistence(data),
                          predictions=plan["predictions"][key], repeated_errors=case["repeated"],
                          reset_errors=case["reset"])
    summary = dict(path=str(path.relative_to(ROOT)), sha256=sha(path),
                   frozen_sha256=record["frozen_prediction_sha256"], cases=cases,
                   holdout_passed=analysis["longer_holdout_passed"], time_step=analysis["time_step"],
                   dry_carbon_mass_kg=plan["dry_carbon_mass_kg"], marker_count=plan["marker_count"],
                   periods=plan["periods"], period_length_m=plan["period_length_m"],
                   carrier_flow_m3_s=plan["carrier_flow_m3_s"])
    write_json(HERE/"reference_summary.json", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    summarize((args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve())
