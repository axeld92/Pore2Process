"""Static scientific comparison of the completed end-to-end runs."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render(path):
    path=Path(path)
    summary=json.loads(path.read_text())
    cases=summary["designs"]
    labels=[f"{c['radius_m']*1e6:.0f} µm radius" for c in cases]
    colors=["#007f86","#cd633f"]
    fig,axes=plt.subplots(1,3,figsize=(11.8,4),layout="constrained")
    values=([c["translation"]["permeability_m2"]*1e12 for c in cases],
            [c["translation"]["pressure_drop_Pa"] for c in cases],
            [c["dwsim"]["shaft_power_W"]*1e3 for c in cases])
    settings=(("OpenLB → Julia","Measured permeability (10⁻¹² m²)"),
              ("Julia Darcy translation","Required pressure rise (Pa)"),
              ("Actual DWSIM pump solve","Shaft power (mW)"))
    for ax,y,(title,ylabel) in zip(axes,values,settings):
        bars=ax.bar(labels,y,color=colors,width=.6)
        ax.bar_label(bars,fmt="%.3g",padding=4)
        ax.set(title=title,ylabel=ylabel,ylim=(0,max(y)*1.2))
        ax.grid(axis="y",alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle("Smaller pores require more pumping power at the same flow\n"
                 "Fixed 200 µm × 200 µm pore cell; idealized 1 m slab, 1 m² area, 0.6 L/min water, 75% pump efficiency",fontsize=11)
    target=path.parent / "design_comparison.png"
    fig.savefig(target,dpi=180)
    plt.close(fig)
    return target


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path,nargs="?")
    path=parser.parse_args().summary
    if path is None:
        from pore2process_demo.paths import ROOT
        root = ROOT
        path=Path(json.loads((root / "results/coupling/latest.json").read_text())["summary"])
    print(render(path))
