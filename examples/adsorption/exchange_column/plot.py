"""Plot the frozen predictions and the two source marker series."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.alhamed.study import observations
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import CASES, HERE, OUTPUT


def render(path, record, frozen):
    base = path.parent
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0), layout="constrained", sharey=True)
    compact = {}
    for ax, (size, (series, diameter)) in zip(axes, CASES.items(), strict=True):
        t, fraction = observations(series)
        ax.scatter(t/60, fraction, facecolors="none", edgecolors="#202020", s=20, label="Figure 7 markers", zorder=4)
        case = frozen["cases"][size]
        compact[size] = {}
        for key, color, style, label in (("primary","#237b78","-","Path exchange"),
                                         ("control","#a05b39","--","Same paths, no exchange")):
            saved = read_json(base/case[key])
            data, inputs = saved["result"], saved["inputs"]
            ax.plot(np.array(data["time_s"])/60, np.array(data["outlet_concentration_kg_m3"])/inputs["inlet_concentration_kg_m3"],
                    color=color, ls=style, lw=1.8, label=label)
            compact[size][key] = dict(inputs=inputs, trace=saved["trace"],
                                      normalized_mass_residual=data["normalized_mass_residual"],
                                      normalized_interface_flux_mismatch=data["normalized_interface_flux_mismatch"])
        ax.set(xlim=(0,300), ylim=(0,1.02), title=f"{1000*diameter:g} mm carbon particles",
               xlabel="Time (min)")
        ax.grid(alpha=.16)
    axes[0].set_ylabel("Outlet / inlet phenol concentration")
    axes[1].legend(fontsize=8, loc="lower right")
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"figure7_exchange.png",dpi=180)
    plt.close(fig)
    write_json(HERE/"reference_summary.json",dict(summary=str(path.relative_to(ROOT)),sha256=sha(path),
        qualified=record["qualified"], figure7=record["figure7"], cases=compact,
        closure_resolution=read_json(base/"closure.json")["attempts"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve()
    record = validate_record(path)
    render(path, record, read_json(path.parent/"frozen_predictions.json"))


if __name__=="__main__":
    main()
