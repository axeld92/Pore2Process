"""Plot the eventual escape and its approach in each synthetic column."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import read_json
from .model import compose, during_holdout
from .run import FLOWS, HERE, STEPS, TARGET, check, inputs, probe_field


def plot(path):
    check(path)
    record = read_json(path)
    _, _, probes = inputs()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), sharey=True, layout="constrained")
    with np.load(path.parent / "clean_film.npz") as film:
        for ax, flow, rate in zip(axes, FLOWS, (23.3, 62.5), strict=True):
            wall = record["perfect_sink"][flow]
            step = STEPS["fine"]
            sections = during_holdout(probe_field(probes, flow, "holdout"), step)
            perfect = compose(sections, wall["blocks"], wall["terminal_sections"])
            dt = step * wall["physical_time_scale"] / wall["liquid_residence_s"]
            cases = (
                (film[f"{flow}_wall"], "Clean film, wall-hit allocation", "#9b5431", "-"),
                (film[f"{flow}_uniform_per_flow"], "Clean film, best allocation", "#5c5285", "--"),
                (perfect, "Perfect first-contact sink", "#176e67", "-."),
            )
            for pulse, label, color, style in cases:
                ax.plot(dt * np.arange(len(pulse)), 100 * np.cumsum(pulse),
                        label=label, color=color, ls=style, lw=1.8)
            ax.axhline(100 * TARGET, color="#777c82", ls=":", label="5% cutoff")
            ax.set(title=f"{rate:g} mL/min", xlim=(0, 2), ylim=(0, 13),
                   xlabel="Time / nominal liquid residence time")
            ax.grid(alpha=.15)
    axes[0].set_ylabel("Inlet dose escaped (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, frameon=False, loc="outside lower center", fontsize=8.8)
    dest = HERE / "figures" / "early_escape.png"
    dest.parent.mkdir(exist_ok=True)
    fig.savefig(dest, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    plot(parser.parse_args().summary.resolve())
