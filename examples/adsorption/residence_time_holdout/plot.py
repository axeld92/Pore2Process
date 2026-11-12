"""Plot the complete held-out rise and tail from the saved cohort."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.porous_outlet_tracer.run import history
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, TAIL_S


def render(path, record):
    base = path.parent
    times, model = np.loadtxt(base / "prediction.csv", delimiter=",", skiprows=1, unpack=True)
    data = history(base / "holdout/history.csv")
    pore = data["recovered_kg"] / record["request"]["dose_kg"]
    measured = record["analysis"]["measured"]
    fig, (rise, tail) = plt.subplots(1, 2, figsize=(9.6, 3.9), layout="constrained")
    for cumulative, label, color, style in ((pore, "Pore simulation", "#237b78", "-"),
                                            (model, "Frozen path model", "#b05a31", "--")):
        rise.plot(1000 * times, cumulative, label=label, color=color, ls=style)
        remaining = 1 - cumulative
        tail.semilogy(times, np.where(remaining > 1e-6, remaining, np.nan),
                      label=label, color=color, ls=style)
    rise.set(title="Arrival and main rise", xlim=(100, 1000 * max(.35, measured["t90"])),
             xlabel="Time from pulse onset (ms)", ylabel="Fraction of injected tracer recovered",
             ylim=(0, 1.02))
    rise.legend(fontsize=9)
    tail.set(title="Late tracer recovery", xlim=(.25, times[-1]), ylim=(1e-6, .3),
             xlabel="Time from pulse onset (s)", ylabel="Fraction of injected tracer still in bed")
    tail.axvline(TAIL_S, color="#555555", lw=.8)
    folder = HERE / "figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder / "outlet_holdout.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT / read_json(OUTPUT / "latest.json")["summary"]).resolve()
    record = validate_record(path)
    render(path, record)
    write_json(HERE / "reference_summary.json",
               dict(summary=str(path.relative_to(ROOT)), summary_sha256=sha(path),
                    model=read_json(path.parent / "frozen_prediction.json")["model"],
                    analysis=record["analysis"]))


if __name__ == "__main__":
    main()
