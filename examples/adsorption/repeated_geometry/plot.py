"""Plot the measured bypass and the frozen eight-cell scenarios."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT


def plot(path):
    recorded = validate_record(path)
    plan = read_json(path.parent/"frozen_prediction.json")
    results = recorded["analysis"]["cases"]
    cells = np.arange(1, plan["periods"]+1)
    fig, axes = plt.subplots(2, 2, sharey="row", sharex=True, figsize=(10, 6.25), layout="constrained")
    for j, (size, diameter) in enumerate((("d800", .80), ("d450", .45))):
        primary = size+"_primary"
        half = size+"_half_step"
        for i, (key, ylabel) in enumerate((("survival", "All tracer wall-free (%)"),
                                          ("fast_survival", "Initially fast tracer wall-free (%)"))):
            ax = axes[i, j]
            held = results[primary]["observed"]
            fine = results[half]["observed"]
            data = plan["predictions"][primary]
            ax.plot(cells, 100*np.array(data["repeated"][key]), color="#237675", lw=2.,
                    label="Frozen repeated model")
            ax.plot(cells, 100*np.array(data["reset_after_two_cells"][key]), color="#a75832",
                    ls="--", lw=1.7, label="Within-class reset every two cells")
            ax.plot(cells, 100*np.array(held[key]), ls="none", marker="o", color="#1d2335",
                    mfc="white", ms=5., label="Tracked at 80 µs")
            ax.plot(cells, 100*np.array(fine[key]), ls="none", marker="x", color="#74649e",
                    ms=5., label="Tracked at 40 µs")
            ax.axvspan(.7, 2.05, color="#cccccc", alpha=.18)
            ax.grid(alpha=.16)
            ax.set(xticks=cells, xlim=(.75, 8.25))
            if j == 0:
                ax.set_ylabel(ylabel)
        axes[0, j].set_title(f"{diameter:.2f} mm particle-matched Péclet number")
        axes[1, j].set_xlabel("Number of porous cells traversed")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, fontsize=9, frameon=False)
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"eight_cell_bypass.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    plot((args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve())
