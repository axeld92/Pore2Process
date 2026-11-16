"""Plot the two frozen repeat-join scenarios and their concentration difference."""
import argparse
from pathlib import Path

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
    if sha(frozen_path)!=record["frozen_prediction_sha256"]:
        raise ValueError("Frozen prediction manifest changed")
    frozen = read_json(frozen_path)
    markers = points()
    fig,axes = plt.subplots(2,2,figsize=(9.2,5.4),sharex=True,layout="constrained")
    for col,(size,(series,diameter)) in enumerate(CASES.items()):
        ax,delta = axes[:,col]
        curves = {}
        for name,label,color,line in (("mixed","Fully mixed joins","#a75033","--"),
                                      ("preserved","Path-preserving joins","#116f6b","-")):
            curve = frozen["cases"][size][name]["selected"]
            if sha(path.parent/curve)!=frozen["curve_sha256"][curve]:
                raise ValueError("Frozen curve changed")
            run = read_json(path.parent/curve)
            t = np.array(run["result"]["time_s"])/60
            c = np.array(run["result"]["outlet_concentration_kg_m3"])/run["inputs"]["inlet_concentration_kg_m3"]
            ax.plot(t,c,label=label,color=color,linestyle=line,lw=1.8)
            curves[name] = (t,c)
        data = markers[series]
        ax.scatter(data["time_s"]/60,data["fraction"],s=18,facecolors="white",edgecolors="#252934",
                   linewidths=.8,label="Alhamed Figure 7",zorder=3)
        t,c = curves["preserved"]
        b,a = curves["mixed"]
        delta.plot(t,100*(c-np.interp(t,b,a)),color="#116f6b",lw=1.7)
        delta.axhline(0,color="#777",lw=.8)
        delta.set_xlabel("Time (min)")
        ax.set(title=f"{1000*diameter:.2f} mm",ylim=(-.02,1.02),xlim=(0,300))
        for axis in (ax,delta):
            axis.grid(alpha=.17)
    axes[0,0].set_ylabel("Outlet / feed concentration")
    axes[1,0].set_ylabel("Join effect (percentage points)")
    handles,labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=3,frameon=False,fontsize=9)
    dest = HERE/"figures/join_comparison.png"
    dest.parent.mkdir(exist_ok=True)
    fig.savefig(dest,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path,nargs="?")
    args = p.parse_args()
    plot(args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"])
