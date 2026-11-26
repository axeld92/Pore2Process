"""Draw the measured peak times and pulse widths from the Julia evidence screen."""
import argparse
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
COLORS = {"FC": "#007570", "CF": "#a5491a"}


def plot(reference, output):
    with reference.open("rb") as stream:
        data = tomllib.load(stream)
    fig, (timing, width) = plt.subplots(1, 2, figsize=(10.4, 4.5))
    for index, rate in enumerate(data["rates"]):
        q = rate["flow_ml_min"]
        for direction, offset, marker in [("FC", -0.19, "o"), ("CF", 0.19, "s")]:
            rows = [r for r in data["curves"] if (r["flow_ml_min"], r["direction"]) == (q, direction)]
            for repeat, row in enumerate(rows):
                x = index + offset + 0.07 * (repeat - (len(rows) - 1)/2)
                peak = sum(row["peak_midrange_time_min"]) / 2
                color = COLORS[direction]
                timing.vlines(x, *row["possible_peak_time_min"], color=color, linewidth=2)
                timing.plot(x, peak, marker, color=color, markersize=6)
                bounds = [q * value for value in row["fwhm_reading_bounds_min"]]
                width.vlines(x, *bounds, color=color, linewidth=2)
                width.plot(x, q * row["fwhm_min"], marker, color=color, markersize=6)
        for axis, value in [(timing, rate["baseline_peak_min"]), (width, q * rate["baseline_fwhm_min"])]:
            axis.plot(index, value, "D", markerfacecolor="white", markeredgecolor="#333333", markersize=6)
    timing.set(title="Peak arrival on the published clock", ylabel="Time (min)", ylim=(0, 460))
    width.set(title="Liquid delivered across the pulse width", ylabel="Flow × width at half height (mL)", ylim=(20, 66))
    for axis in (timing, width):
        axis.set(xlabel="Flow rate (mL/min)", xlim=(-0.48, 2.48), xticks=range(3), xticklabels=["0.3", "0.4", "1.0"])
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#dce2e1", linewidth=0.65)
        axis.set_axisbelow(True)
    legend = [Line2D([], [], marker=m, linestyle="none", color=c, label=l, markerfacecolor=f)
              for m, c, l, f in [("o", COLORS["FC"], "Fine then coarse", COLORS["FC"]),
                                  ("s", COLORS["CF"], "Coarse then fine", COLORS["CF"]),
                                  ("D", "#333333", "Two-layer calculation", "white")]]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.008))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.24, top=0.9, wspace=0.27)
    fig.savefig(output, dpi=175, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=HERE / "reference.toml")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.reference, args.output)
