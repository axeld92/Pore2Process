"""Plot the outlet comparison and the primary path-exchange matrix."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.porous_outlet_tracer.run import PULSE_S, history
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import ARCHIVE, HERE, OUTPUT


def render(path, record, cases):
    data = history(ARCHIVE / "holdout/history.csv")
    prediction = history(path.parent / "prediction.csv")
    time = prediction["time_s"]
    actual = data["recovered_kg"] / (record["mesh"]["flow_m3_s"] * PULSE_S)
    fig, (rise, tail, matrix) = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
    curves = ((actual, "Pore reference", "#237b78", "-"),
              (prediction["primary"], "With exchange", "#694499", "--"),
              (prediction["persistent"], "Same path each period", "#b05a31", ":"))
    for cdf, label, color, style in curves:
        rise.plot(1000 * time, cdf, label=label, color=color, ls=style)
        remaining = 1 - cdf
        tail.semilogy(time, np.where(remaining > 1e-6, remaining, np.nan), color=color, ls=style)
    rise.set(title="Arrival and main rise", xlim=(125, 350), ylim=(0, 1.01),
             xlabel="Time from pulse onset (ms)", ylabel="Fraction recovered")
    rise.legend(fontsize=8)
    tail.set(title="Late recovery", xlim=(.25, 1.65), ylim=(1e-5, .3),
             xlabel="Time from pulse onset (s)", ylabel="Fraction still in bed")
    counts = np.asarray(cases["primary"]["decile_counts"])
    probability = counts / counts.sum(axis=1, keepdims=True)
    heat = matrix.imshow(probability, cmap="Purples", vmin=0, vmax=1, interpolation="nearest")
    ticks = [0, 4, 9]
    matrix.set(title="Path exchange from markers", xlabel="Second period transit", ylabel="First period transit",
               xticks=ticks, yticks=ticks, xticklabels=["Fastest", "5th", "Slowest"],
               yticklabels=["Fastest", "5th", "Slowest"])
    fig.colorbar(heat, ax=matrix, fraction=.046, label="Conditional probability")
    folder = HERE / "figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder / "path_exchange.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT / read_json(OUTPUT / "latest.json")["summary"]).resolve()
    record = validate_record(path)
    frozen = read_json(path.parent / "frozen_prediction.json")
    render(path, record, frozen["cases"])
    write_json(HERE / "reference_summary.json",
               dict(summary=str(path.relative_to(ROOT)), sha256=sha(path),
                    cases=frozen["cases"], analysis=record["analysis"]))


if __name__ == "__main__":
    main()
