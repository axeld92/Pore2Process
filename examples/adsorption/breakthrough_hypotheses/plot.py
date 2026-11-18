"""Plot the frozen 25 cm hypotheses and corresponding published markers."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.carbon_access_transfer import source
from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT

from .screen import HERE


def run_directory(reference):
    file = ROOT / reference["path"]
    if sha(file) != reference["sha256"]:
        raise ValueError("Expected the frozen results named in the compact record")
    return file.parent


def plot(output):
    record = read_json(HERE / "reference_summary.json")["runs"]
    periodic = run_directory(record["periodic"])
    periodic_cases = read_json(periodic / "frozen.json")["cases"]
    equilibria = run_directory(record["freundlich_column"])
    joint = run_directory(record["joint"]) if "joint" in record else None
    joint_cases = read_json(joint / "frozen.json")["cases"] if joint else {}
    data = source.markers()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.1), sharey=True, layout="constrained")
    for axis, rate, flow in zip(axes, (23.3, 62.5), ("q233", "q625"), strict=True):
        case = f"{flow}_l25"
        experiment = data[f"f4_{case}"]
        confined = run_directory(record[f"confined_{flow}"])
        confined_cases = read_json(confined / "frozen.json")["cases"]
        traces = [
            (periodic / periodic_cases[case]["baseline"]["selected"], "Periodic, Langmuir", "#545962", "-"),
            (confined / confined_cases[case]["wall_flow"]["selected"], "Confined wall-flow scenario", "#286a86", "-"),
            (confined / confined_cases[case]["wall_and_histories"]["selected"], "Confined, separate carbon histories", "#286a86", "--"),
            (equilibria / read_json(equilibria / f"{flow}_levels.json")["selected"], "Periodic, Freundlich", "#a75226", "-"),
        ]
        if joint:
            traces.append((joint / joint_cases[case]["joint"]["selected"], "Confined, adjusted equilibrium", "#287044", "-."))
        for file, label, color, line in traces:
            curve = read_json(file)
            time = np.array(curve["result"]["time_s"]) / 60
            sample = np.linspace(0., experiment["time"][-1] / 60, 1401)
            fraction = np.array(curve["result"]["outlet_concentration_kg_m3"]) / curve["inputs"]["inlet_concentration_kg_m3"]
            axis.plot(sample, 100 * np.interp(sample, time, fraction), label=label, color=color, ls=line, lw=1.6)
        axis.plot(experiment["time"] / 60, 100 * experiment["fraction"], "o", label="Published markers",
                  mfc="white", mec="#202124", ms=4.0)
        axis.set(title=f"25 cm bed, {rate:g} mL/min", xlabel="Time (min)", ylim=(0, 102))
        axis.grid(alpha=.20)
    axes[0].set_ylabel("Outlet / feed concentration (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=3, fontsize=8.7, frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "figures" / "figure4_hypotheses.png")
    plot(parser.parse_args().output)
