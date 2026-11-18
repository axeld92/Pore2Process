"""Show early leakage and the full frozen Figure 7 curves."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.alhamed.source import points
from pore2process_demo.io import read_json, sha

from .run import CASES, HERE

STYLES = {
    "gac4": ("GAC4 Langmuir", "#286a86", "--"),
    "langmuir": ("Alhamed Langmuir", "#545962", "-"),
    "blend": ("50/50 equilibrium blend", "#287044", "-."),
    "freundlich": ("Alhamed Freundlich", "#a75226", "-"),
}


def plot(base, output):
    frozen = read_json(base / "frozen.json")
    if sha(base / "frozen.json") != read_json(base / "summary.json")["frozen_sha256"]:
        raise ValueError("Expected a sealed Figure 7 forecast")
    data = points()
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.5), layout="constrained")
    for (early, full), (size, (series, diameter)) in zip(axes, CASES.items(), strict=True):
        observed = data[series]
        early_peak = 5.
        for name, (label, color, line) in STYLES.items():
            filename = frozen["cases"][size][name]["selected"]
            file = base / filename
            if sha(file) != frozen["curve_sha256"][filename]:
                raise ValueError("Frozen plotted curve changed")
            saved = read_json(file)
            time = np.asarray(saved["result"]["time_s"]) / 60
            feed = saved["inputs"]["inlet_concentration_kg_m3"]
            fraction = 100 * np.asarray(saved["result"]["outlet_concentration_kg_m3"]) / feed
            early_peak = max(early_peak, np.max(fraction[time <= 30]))
            for axis in (early, full):
                axis.plot(time, fraction, label=label, color=color, ls=line, lw=1.5)
        for axis in (early, full):
            axis.plot(observed["time_s"] / 60, 100 * observed["fraction"], "o", label="Figure 7 markers",
                      ms=4, mfc="white", mec="#202124", zorder=5)
            axis.grid(alpha=.20)
            axis.set_xlabel("Time (min)")
        early.axhline(5., color="#535862", lw=.8, ls=":")
        early.set(xlim=(0, 30), ylim=(0, min(101, 10*np.ceil((early_peak+1)/10))), ylabel="Outlet / feed (%)",
                  title=f"{diameter * 1000:.2f} mm, early breakthrough")
        full.set(xlim=(0, observed["time_s"][-1]/60 * 1.02), ylim=(0, 101),
                 title=f"{diameter * 1000:.2f} mm, full curve")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=3, fontsize=8.8, frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("--output", type=Path, default=HERE/"figures/figure7_equilibrium.png")
    args = parser.parse_args()
    plot(args.base, args.output)
