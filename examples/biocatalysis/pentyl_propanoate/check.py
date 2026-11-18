"""Check retained evidence and replay the held-out comparison without simulators."""

import json

from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT

from .run import HERE, compare, read_rows


def check():
    record = read_json(HERE / "reference/summary.json")
    bed = record["bed"]
    inputs = record["batch"]["sources"] | {
        record["held_out_conversion_source"]: record["held_out_conversion_sha256"],
        bed["mixing_source"]: bed["mixing_source_sha256"],
        str((HERE / "model.jl").relative_to(ROOT)): record["model_sha256"],
        str((HERE / "run.py").relative_to(ROOT)): record["runner_sha256"],
    }
    for name, digest in inputs.items():
        if sha(ROOT / name) != digest:
            raise ValueError(f"Changed source for retained comparison: {name}")
    table = HERE / "reference/bed_comparison.csv"
    if sha(table) != record["bed_comparison_sha256"]:
        raise ValueError("Retained held-out comparison table changed.")
    rows = read_rows(table)
    updated, metrics = compare(rows)
    if any(float(a["observed_acid_conversion"]) != b["observed_acid_conversion"]
           for a, b in zip(rows, updated, strict=True)):
        raise ValueError("A retained observation differs from the source extraction.")
    for model, values in metrics.items():
        for metric, value in values.items():
            if abs(value - bed["comparison"][model][metric]) > 1e-10:
                raise ValueError(f"Retained metric differs: {model}/{metric}")
    return metrics


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
