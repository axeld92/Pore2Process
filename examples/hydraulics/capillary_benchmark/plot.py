"""Render scientific comparisons from a completed hydraulic result snapshot."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render(summary_path):
    summary_path = Path(summary_path).resolve()
    summary = json.loads(summary_path.read_text())
    case = summary["cases"][2]
    directory = summary_path.parent / case["name"]
    raw = case["inputs_and_measurements"]
    planes = np.genfromtxt(directory / "planes.csv", delimiter=",", names=True)
    sample = np.genfromtxt(directory / "midplane.csv", delimiter=",", names=True)
    radius = raw["diameter_m"]/2
    r = np.hypot(sample["y_m"]-radius, sample["z_m"]-radius)
    radial = np.linspace(0,radius,150)
    analytical_u = raw["imposed_delta_p_Pa"]/(4*raw["viscosity_Pa_s"]*raw["length_m"])*(radius**2-radial**2)
    fig, axes = plt.subplots(1,3,figsize=(13.3,3.8),layout="constrained")
    axes[0].plot(radial*1e6,analytical_u*1e3,label="Hagen–Poiseuille",color="#444444")
    axes[0].scatter(r*1e6,sample["ux_m_s"]*1e3,s=9,label=f"OpenLB, N={raw['resolution']}",alpha=.5,color="#006d77")
    axes[0].set(xlabel="Distance from pore axis (µm)",ylabel="Axial velocity (mm/s)",title="Midplane velocity")
    axes[0].legend(fontsize=8)
    axes[1].plot(planes["x_m"]*1e3,planes["mean_pressure_Pa"],label="OpenLB",color="#006d77")
    axes[1].plot(planes["x_m"]*1e3,raw["imposed_delta_p_Pa"]*(1-planes["x_m"]/raw["length_m"]),"--",label="Analytical",color="#444444")
    axes[1].set(xlabel="Axial position x (mm)",ylabel="Mean gauge pressure (Pa)",title="Pressure profile")
    axes[1].legend(fontsize=8)
    grid=sorted((c for c in summary["cases"] if c["name"] in ("n16","n24","n32","n40")),
                key=lambda c:c["inputs_and_measurements"]["resolution"])
    n=[c["inputs_and_measurements"]["resolution"] for c in grid]
    error=[100*c["metrics"]["permeability_relative_error"] for c in grid]
    axes[2].plot(n,error,"o-",color="#006d77")
    axes[2].axhline(1,color="#999999",linestyle="--",label="1% acceptance")
    axes[2].set(xlabel="Cells across pore diameter",ylabel="Permeability error (%)",title="Grid refinement",xticks=n)
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    target=summary_path.parent / "hydraulic_benchmark.png"
    fig.savefig(target,dpi=180)
    plt.close(fig)
    return target


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path,nargs="?")
    selected=parser.parse_args().summary
    if selected is None:
        from pore2process_demo.paths import ROOT
        root = ROOT
        selected=Path(json.loads((root / "results/hydraulic/latest.json").read_text())["summary"])
    print(render(selected))
