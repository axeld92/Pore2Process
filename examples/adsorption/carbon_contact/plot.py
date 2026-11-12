"""Plot modeled carbon allocation and observed wall contact by transit rank."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .prepare import HERE, OUTPUT

PRIMARY = (
    ("reference_primary", "Saved-flow reference", "#44759a"),
    ("d800_primary", "0.80 mm Péclet match", "#176b62"),
    ("d450_primary", "0.45 mm Péclet match", "#785596"),
)
PREVIOUS = (
    ("void", "Carbon per void volume", "#a6abb0"),
    ("slow", "Strong slow-path enrichment", "#bb653b"),
)
WALL = (
    ("Only slowest group observed", "#705197"),
    ("Slowest 10% and other ranks", "#447b9b"),
    ("Other ranks only", "#8cb2a5"),
    ("No recorded hit", "#e6e8ea"),
)


def _mass_bars(ax, data):
    scenarios = data["previous_carbon_scenarios"]
    cases = data["cases"]
    bars = [
        (label, color, [scenarios[name]["fast_fraction"], scenarios[name]["slow_fraction"]])
        for name, label, color in PREVIOUS
    ]
    for name, label, color in PRIMARY:
        fractions = np.asarray(cases[name]["advection_ranked_exposure"]["patch_carbon_fraction"])
        bars.append((label, color, fractions[[0, -1]]))

    centres = np.arange(2)
    width = .156
    for j, (label, color, share) in enumerate(bars):
        x = centres + width * (j - (len(bars) - 1) / 2)
        ax.bar(x, 100 * np.asarray(share), width=width * .94, color=color, label=label, zorder=2)

    ax.set(title="A   Contact-based carbon allocation", ylabel="Share of fixed dry carbon (proxy)",
           xlabel="Transit-time group (advection rank for marker cases)",
           xticks=centres, xticklabels=["Fastest 10%", "Slowest 10%"],
           xlim=(-.5, 1.5), ylim=(0, max(54, 110 * max(max(share) for _, _, share in bars))))
    ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.grid(axis="y", color="#cbd1d4", linewidth=.65)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, fontsize=8, handlelength=1.15)


def _wall_bars(ax, data):
    rows = []
    for name, _, _ in PRIMARY:
        stats = data["pooled_face_overlap"][name.removesuffix("_primary")]
        slow = stats["slow_only_observed_area"]
        shared = stats["shared_slow_and_other_area"]
        visited = stats["area_visited"]
        rows.append([slow, shared, visited - slow - shared, 1 - visited])
    coverage = np.asarray(rows)
    np.testing.assert_allclose(coverage.sum(axis=1), 1., rtol=0, atol=1e-12)
    if np.min(coverage) < -1e-12:
        raise ValueError("Wall coverage categories must be nonnegative")

    y = np.arange(len(rows))
    left = np.zeros(len(rows))
    for values, (label, color) in zip(coverage.T, WALL, strict=True):
        width = 100 * values
        ax.barh(y, width, left=left, height=.43, color=color, label=label,
                edgecolor="white", linewidth=.5, zorder=2)
        left += width
    ax.set(title="B   Observed wall coverage by realized transit rank", xlabel="Share of corrected grain surface area",
           yticks=y, yticklabels=[label for _, label, _ in PRIMARY],
           xlim=(0, 100), ylim=(2.55, -.55))
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.grid(axis="x", color="#cbd1d4", linewidth=.65)
    ax.set_axisbelow(True)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.25), ncol=2,
              frameon=False, fontsize=8, handlelength=1.15)


def render(data, output):
    fig, (mass, walls) = plt.subplots(1, 2, figsize=(11.25, 4.1), layout="constrained")
    _mass_bars(mass, data)
    _wall_bars(walls, data)
    for ax in (mass, walls):
        ax.spines[["top", "right"]].set_visible(False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT / read_json(OUTPUT / "latest.json")["summary"]).resolve()
    record = validate_record(path)
    data = record["analysis"]
    render(data, HERE / "figures" / "carbon_contact.png")
    selected = {name:data["cases"][name] for name,_,_ in PRIMARY}
    write_json(HERE/"reference_summary.json",dict(summary=str(path.relative_to(ROOT)),sha256=sha(path),
        total_dry_mass_kg=data["total_dry_mass_kg"],wall=data["wall"],cases=selected,
        pooled_face_overlap=data["pooled_face_overlap"],step_refinement=data["step_refinement"],
        independent_sample=data["independent_sample"],checks_passed=data["all_contact_checks_passed"]))


if __name__ == "__main__":
    main()
