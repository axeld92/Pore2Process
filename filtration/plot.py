"""Scientific figure for the saved fixed-flux, fixed-selective-layer comparison."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render(path):
    path=Path(path)
    data=json.loads(path.read_text())
    ratio=data["supports"][0]["permeability_m2"]/data["supports"][1]["permeability_m2"]
    alpha=np.geomspace(.01,1000,300)
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout="constrained")
    blue,orange="#007c91","#cf6136"
    ax=axes[0]
    ax.semilogx(alpha,100*(1-ratio)/(1+alpha),color=blue,label="Series-law prediction")
    for s in data["scenarios"]:
        ax.plot(s["alpha"],100*s["shaft_energy_saving_fraction"],"o",color=blue)
    ax.plot([],[],"o",color=blue,label="Actual DWSIM pairs")
    ax.set(xlabel="Selective resistance / baseline support resistance, α",
           ylabel="Shaft energy saving (%)",ylim=(0,65),title="Benefit is limited by the support's share")
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    ax=axes[1]
    ax.semilogx(alpha,100*alpha/(1+alpha),color=blue,label="4 µm support pores")
    ax.semilogx(alpha,100*alpha/(ratio+alpha),color=orange,label="5 µm support pores")
    ax.axhline(50,color="grey",lw=.8,ls="--");ax.axhline(90,color="grey",lw=.8,ls=":")
    ax.set(xlabel="Selective resistance / baseline support resistance, α",
           ylabel="Selective share of total resistance (%)",ylim=(0,102),title="Define dominance: >50%; strong: ≥90%")
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[2]
    x=np.array([s["alpha"] for s in data["scenarios"]])
    for i,color in enumerate((blue,orange)):
        energy=[s["cases"][i]["actual_specific_shaft_energy_kWh_m3"] for s in data["scenarios"]]
        ax.loglog(x,np.array(energy)*1e6,"o-",color=color,label=f"{4+i} µm pores")
    ax.set_xticks(x,[f"α = {s['alpha']:g}" for s in data["scenarios"]])
    ax.set(ylabel="Actual DWSIM shaft energy (10⁻⁶ kWh/m³)",title="Same flux and selective layer within each pair")
    ax.legend(fontsize=8);ax.grid(alpha=.2);ax.set_axisbelow(True)
    settings=data["settings"]
    fig.suptitle(f"Synthetic clean-water support screening: {settings['flux_LMH']:g} L m⁻² h⁻¹, "
        f"{settings['support_thickness_m']*1e3:g} mm support, 75% shaft-to-hydraulic efficiency\n"
        "OpenLB permeability measurements → Julia layers in series → saved DWSIM pump cases; no experimental calibration",fontsize=10)
    target=path.parent / "filtration_comparison.png"
    fig.savefig(target,dpi=180);plt.close(fig)
    return target


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("summary",type=Path)
    print(render(parser.parse_args().summary))
