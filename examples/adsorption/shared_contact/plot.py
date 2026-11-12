"""Plot frozen shared-contact predictions alongside the prior column comparison."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.alhamed.study import observations
from examples.adsorption.contact_column.run import scoring_curve
from examples.adsorption.exchange_column.run import CASES
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE, OUTPUT, PREVIOUS

STYLE = {
    "wall": ("Shared carbon, wall-weighted area","#237b78","-"),
    "flow": ("Shared carbon, flow-weighted area","#694499","-."),
}


def render(path,record):
    base = path.parent
    frozen = read_json(base/"frozen_predictions.json")
    previous = read_json(PREVIOUS/"frozen_predictions.json")
    fig, axes = plt.subplots(2,2,figsize=(10,6.3),layout="constrained",sharey="row",
                             height_ratios=(1.45,1))
    for column,(size,(series,diameter)) in enumerate(CASES.items()):
        a,b = axes[:,column]
        t,y = observations(series)
        a.scatter(t/60,y,facecolors="none",edgecolors="#202020",s=18,label="Figure 7 markers",zorder=4)
        early = t<=60
        b.scatter(t[early],y[early],facecolors="none",edgecolors="#202020",s=18,zorder=4)
        old = previous["cases"][size]["advection_rank"]
        old_result = read_json(PREVIOUS/old["selected"])["result"]
        curves = [(old_result,"Separate carbon, wall-weighted","#9f6037","--")]
        for name,item in frozen["cases"][size].items():
            full = read_json(base/item["full"][-1])["result"]
            short = read_json(base/item["short"][-1])["result"] if item["short"] else None
            curves.append((scoring_curve(full,short),*STYLE[name]))
        for result,label,color,style in curves:
            time = np.array(result["time_s"])
            c = np.array(result["outlet_concentration_kg_m3"])/.15
            a.plot(time/60,c,label=label,color=color,ls=style,lw=1.65)
            b.plot(time,c,color=color,ls=style,lw=1.65)
        a.set(xlim=(0,300),ylim=(0,1.02),title=f"{1000*diameter:g} mm particles",xlabel="Time (min)")
        b.set(xlim=(0,60),ylim=(0,1.02),xlabel="First minute (s)")
        for ax in (a,b):
            ax.grid(alpha=.17)
    for ax in axes[:,0]:
        ax.set_ylabel("Outlet / inlet phenol concentration")
    handles,labels = axes[0,1].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=2,fontsize=8,frameon=False)
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"figure7_shared_contact.png",dpi=180)
    plt.close(fig)
    write_json(HERE/"reference_summary.json",dict(summary=str(path.relative_to(ROOT)),sha256=sha(path),
        figure7=record["figure7"],checks={size:{name:item["review"] for name,item in cases.items()}
                                       for size,cases in frozen["cases"].items()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve()
    render(path,validate_record(path))


if __name__=="__main__":
    main()
