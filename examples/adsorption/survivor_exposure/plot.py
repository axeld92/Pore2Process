"""Plot the held-out survivor capture and the frozen column predictions."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.alhamed.study import observations
from examples.adsorption.exchange_column.run import CASES
from pore2process_demo.io import read_json
from pore2process_demo.paths import ROOT

from .column import curve
from .run import ARCHIVE, HERE, OUTPUT


def render(path):
    record = read_json(path)
    holdout = read_json(ROOT/record["holdout"])["scorecard"]
    frozen = read_json(path.parent/"frozen_predictions.json")
    prior = read_json(ARCHIVE/"frozen_predictions.json")
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.6), layout="constrained", sharey="row")
    for i, (size, (series, diameter)) in enumerate(CASES.items()):
        sink, outlet = axes[:, i]
        primary = holdout["cases"][size+"_primary"]["models"]
        reference = primary["twin_quartiles"]
        observed = [reference["observed_reset_total_capture"],
                    reference["observed_survivor_conditional_capture"]]
        predicted = [reference["predicted_reset_total_capture"],
                     reference["predicted_survivor_conditional_capture"]]
        x = np.arange(2)
        sink.bar(x-.17, 100*np.array(observed), .34, label="Tracked tracer", color="#456d8b")
        sink.bar(x+.17, 100*np.array(predicted), .34, label="Frozen prediction", color="#d78d50")
        sink.set(xticks=x, xticklabels=["Fresh tracer", "First-cell survivors"], ylim=(0, 5.2),
                 title=f"{diameter*1000:g} mm particles")
        t, y = observations(series)
        outlet.scatter(np.asarray(t)/60, y, facecolors="white", edgecolors="#171717", s=24,
                       label="Figure 7", zorder=5)
        for folder, files, label, color, style in (
                (path.parent, frozen["cases"][size]["quartiles"], "Exposure subgroups", "#17726e", "-"),
                (path.parent, frozen["cases"][size]["parent_only"], "Transit classes", "#575286", "--"),
                (ARCHIVE, prior["cases"][size], "Previous wall-uptake films", "#a65b30", ":")):
            data = curve(folder, files)
            outlet.plot(np.array(data["time_s"])[1:]/60,
                        np.array(data["outlet_concentration_kg_m3"])[1:]/.15,
                        label=label, color=color, linestyle=style, lw=1.75)
        outlet.set(xscale="log", xlim=(.045, 310), ylim=(0, 1.01), xlabel="Time (min, log scale)",
                   xticks=[.05, .1, 1, 10, 100, 300], xticklabels=["0.05", "0.1", "1", "10", "100", "300"])
        sink.grid(axis="y", alpha=.18)
        outlet.grid(alpha=.18)
    axes[0, 0].set_ylabel("Second-cell wall capture (%)")
    axes[1, 0].set_ylabel("Outlet / inlet phenol concentration")
    axes[0, 0].legend(fontsize=8, frameon=False, loc="upper right")
    handles, labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8, frameon=False)
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"survivor_exposure_figure7.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    path = args.summary or ROOT/read_json(OUTPUT/"latest_column.json")["summary"]
    render(path.resolve())
