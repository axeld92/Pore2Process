"""Plot frozen breakthrough predictions and the Figure 7 experimental markers."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.alhamed.source import points
from pore2process_demo.io import read_json, sha
from .run import CASES, HERE, OUTPUT, ROOT


def plot(path):
    record = read_json(path)
    frozen_path = path.parent/"frozen_predictions.json"
    if sha(frozen_path)!=record["frozen_predictions_sha256"]:
        raise ValueError("Frozen curves changed")
    frozen = read_json(frozen_path)
    markers = points()
    fig,axes = plt.subplots(2,2,figsize=(9.2,6.0),layout="constrained")
    legend = {}
    styles = (
        ("selected","Travel-time memory and wall contact","#136e68","-"),
        ("flow_area","Area proportional to flow","#9c512f","--"),
        ("one_age_group","One age group per path","#5c5782",":"),
    )
    for col,(size,(series,diameter)) in enumerate(CASES.items()):
        case = frozen["cases"][size]
        data = markers[series]
        for row in (0,1):
            ax = axes[row,col]
            for key,label,color,line in styles:
                name = case["selected"] if key=="selected" else case["controls"][key]
                if sha(path.parent/name)!=frozen["curves_sha256"][name]:
                    raise ValueError(f"Frozen prediction changed: {name}")
                run = read_json(path.parent/name)
                result,p = run["result"],run["inputs"]
                t = np.asarray(result["time_s"])/60
                c = np.asarray(result["outlet_concentration_kg_m3"])/p["inlet_concentration_kg_m3"]
                artist, = ax.plot(t,c,color=color,ls=line,lw=1.8,label=label)
                legend[label] = artist
            label = "Alhamed Figure 7"
            legend[label] = ax.scatter(data["time_s"]/60,data["fraction"],s=23,facecolors="white",
                edgecolors="#252934",linewidths=0.9,zorder=3,label=label)
            ax.grid(alpha=.18)
            ax.set_ylim(-.025,1.025)
        axes[0,col].set(title=f"{diameter*1e3:.2f} mm",xlim=(0,300),xlabel="Time (min)")
        axes[1,col].set(xscale="log",xlim=(.02,50),xlabel="Time (min, log scale)")
        axes[1,col].axhline(.05,color="#888",lw=.8)
    fig.supylabel("Outlet concentration / feed concentration")
    fig.legend(legend.values(),legend.keys(),loc="outside lower center",ncol=2,frameon=False,fontsize=8.8)
    dest = HERE/"figures/figure7_memory_column.png"
    dest.parent.mkdir(exist_ok=True)
    fig.savefig(dest,dpi=180)
    plt.close(fig)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path,nargs="?")
    args = p.parse_args()
    plot(args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"])
