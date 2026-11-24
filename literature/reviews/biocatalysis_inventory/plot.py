"""Plot the audited source readings; no raw papers or reactor holdouts needed."""
import argparse
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCES = HERE.parents[1] / "sources"


def read(path):
    with path.open("rb") as stream:
        return tomllib.load(stream)


def axis(px, p0, p1, a0, a1):
    return a0 + (np.asarray(px) - p0) * (a1 - a0) / (p1 - p0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=HERE / "reference.toml")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = read(args.summary)
    colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#E69F00"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    fig.suptitle("Independent catalyst inventory and activity evidence", fontsize=15)

    a = axs[0, 0]
    f = read(SOURCES / "damstrup2007" / "evidence.toml")["figure2"]
    for s, fill, color, marker in zip(
        f["series"], summary["damstrup"]["wetting"], colors, ["o", "s", "^"]
    ):
        t = axis(s["x_px"], f["x0_px"], f["x40_px"], 0, 40)
        x = axis(s["y_px"], f["y0_px"], f["y60_px"], 0, 60)
        label = f"{s['mass_g']:g} g, free-wet layer / column {100*fill['free_wet_layer_fraction']:.0f}%"
        a.plot(t, x, "-"+marker, color=color, label=label, markersize=4)
    a.set(title="Damstrup: unequal wet filling at matched W/Q",
          xlabel="Nominal catalyst-based time (min)", ylabel="MAG content (wt%)")
    a.legend(frameon=False, fontsize=8.5)

    a = axs[0, 1]
    t_lin = np.linspace(0, 30, 300)
    for s, color in zip(summary["rejasse"]["storage_series"], ["#0072B2", "#D55E00"]):
        t = np.asarray(s["time_min"])
        y = np.asarray(s["observed_activity_percent"])
        a.scatter(t[:-1], y[:-1], s=25, color=color, label=s["substrate"])
        a.scatter(t[-1:], y[-1:], s=60, facecolor="white", edgecolor=color,
                  linewidth=1.8, zorder=5)
        for model in s["shape_fits"]:
            f = model["retained_fraction"]
            k = model["k_natural_per_min"]
            a.plot(t_lin, 100*(f+(1-f)*np.exp(-k*t_lin)), color=color,
                   linestyle="--" if f == 0 else "-", alpha=0.85)
    a.axvline(10, color="0.65", lw=1)
    a.text(10.8, 98, "30 min open markers held out\nSolid: retained fraction\nDashed: single exponential",
           fontsize=8.5, va="top")
    a.set(title="Réjasse: substrate exposure at 100 °C",
          xlabel="Preincubation time (min)", ylabel="Residual activity (% fresh)", ylim=(0, 107))
    a.legend(frameon=False, fontsize=8.5, loc="center right", bbox_to_anchor=(1, 0.55))

    a = axs[1, 0]
    f = summary["lerin"][2]
    a.plot(f["cycle"], f["activity_percent"], "s-", color=colors[0], label="Separate activity assay")
    a.plot(f["cycle"], f["conversion_percent"], "o-", color=colors[2], label="Reaction conversion")
    a.set(title="Lerin: glyceryl benzoate reuse",
          xlabel="Cycle number (not continuous time)", ylabel="Reported percentage", xticks=range(1, 11))
    a.legend(frameon=False, fontsize=9)

    a = axs[1, 1]
    f = read(SOURCES / "strnisa2018" / "evidence.toml")["figure5"]
    for color, s in zip(colors, f["series"]):
        x = axis(s["y_px"], s["y0_px"], s["y100_px"], 0, 100)
        a.plot(s["tau_min"], x, "o-", color=color, markersize=3.5,
               lw=2 if s["reactor"] == "A" else 1, label=s["reactor"])
    a.set(title="Strniša: scale-up and early-conversion differences",
          xlabel="Calculated liquid residence time (min)", ylabel="Butanol conversion (%)")
    a.legend(ncols=3, frameon=False, fontsize=8.5, title="Reactor", loc="lower right")

    for a in axs.flat:
        a.grid(alpha=0.18)
        a.spines[["top", "right"]].set_visible(False)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
