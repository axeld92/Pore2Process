"""Plot frozen breakthrough forecasts and the observed Figure 4 markers."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.carbon_access_transfer import source
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json

from .run import BASELINE, FLOWS, HERE, PRIMARY, archive_curve


def curve(saved):
    result = saved["result"]
    return (np.asarray(result["time_s"]) / 60,
            np.asarray(result["outlet_concentration_kg_m3"]) /
            saved["inputs"]["inlet_concentration_kg_m3"] * 100)


def plot(base):
    record = validate_record(base / "summary.json")
    plan = read_json(base / "plan.json")
    frozen = read_json(base / "frozen.json")
    markers = source.markers()
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 5.8), layout="constrained")
    styles = (
        (BASELINE, "Previous film correlation", "#737982", ":"),
        ("film45_path45", "N45 film and transport", "#a05b34", "--"),
        (PRIMARY, "N60 film and transport", "#176e67", "-"),
    )
    for (early, full), flow in zip(axes, FLOWS, strict=True):
        condition = plan["cases"][flow]
        key = f"f4_{flow}_l25"
        data = markers[key]
        comparison = record["score"]["comparisons"][key]
        profiles = {name: curve(read_json(base / item["selected"]))
                    for name, item in frozen["cases"][flow].items()}
        profiles[BASELINE] = curve(archive_curve(condition["baseline"]))
        time = np.linspace(0., data["time"][-1] / 60, 2501)
        family = [np.interp(time, t, y) for key, (t, y) in profiles.items() if key != BASELINE]
        crossings = [item["t05_s"]["predicted_s"] for item in comparison["variants"].values()]
        limit = 1.10 * max(x for x in [comparison["observed"]["t05_s"], *crossings, 60.] if x is not None) / 60
        for ax in (early, full):
            ax.fill_between(time, np.min(family, axis=0), np.max(family, axis=0),
                            color="#176e67", alpha=.12, label="Range across the four combinations")
            for name, label, color, style in styles:
                ax.plot(*profiles[name], label=label, color=color, ls=style, lw=1.65)
            ax.plot(data["time"] / 60, data["fraction"] * 100, "o", label="Figure 4 markers",
                    ms=3.5, mfc="white", mec="#34393f", mew=.85, zorder=5)
            ax.grid(alpha=.16)
        early.axhline(5., ls=":", lw=1., color="#444850")
        early.set(xlim=(0, limit), ylim=(0, 14), ylabel="Outlet / feed (%)",
                  title=f'{condition["flow_mL_min"]:g} mL/min · first 5%')
        full.set(xlim=(0, time[-1] * 1.01), ylim=(0, 101),
                 title=f'{condition["flow_mL_min"]:g} mL/min · full curve')
    for ax in axes[-1]:
        ax.set_xlabel("Time (min)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, fontsize=8.6, frameon=False, loc="outside lower center")
    plot_path = HERE / "figures" / "figure4_breakthrough.png"
    plot_path.parent.mkdir(exist_ok=True)
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    plot(parser.parse_args().run.resolve())
