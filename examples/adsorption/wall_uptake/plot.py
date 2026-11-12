"""Show the absorbing-wall check and resulting Figure 7 predictions."""
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

from .estimate import SHARED
from .prepare import HERE, OUTPUT


def render(path,record):
    base = path.parent
    model = read_json(base/"model.json")
    frozen = read_json(base/"frozen_predictions.json")
    prior = read_json(SHARED/"frozen_predictions.json")
    w = np.array(model["closure"]["flow_fractions"])
    rank = 100*(np.cumsum(w)-w/2)
    fig, axes = plt.subplots(2,2,figsize=(10.3,6.5),layout="constrained",sharey="row")
    for column,(size,(series,diameter)) in enumerate(CASES.items()):
        sink,outlet = axes[:,column]
        a = model["sizes"][size]["uptake"]["primary"]
        sink.plot(rank,100*np.array(a["first"]["capture"]),"o",ms=3.2,color="#237b78",label="First-cell capture")
        sink.plot(rank,100*np.array(a["second"]["capture"]),"x",ms=4,color="#524775",label="Second-cell capture")
        sink.plot(rank,-100*np.expm1(-np.array(a["first"]["NTU"])),ls="--",color="#a05e36",
                  label="Constant-conductance prediction")
        sink.set(xlim=(55,100.4),ylim=(-2,102),title=f"{1000*diameter:g} mm particles",
                 xlabel="Advection transit-time percentile")
        t,y = observations(series)
        outlet.scatter(t/60,y,facecolors="none",edgecolors="#202020",s=19,label="Figure 7 markers",zorder=4)
        curves = []
        for archive,item,name,color,ls in (
                (base,frozen["cases"][size],"Wall-uptake films","#237b78","-"),
                (SHARED,prior["cases"][size]["wall"],"Correlation films, wall-weighted","#a05e36","--"),
                (SHARED,prior["cases"][size]["flow"],"Correlation films, flow-weighted","#524775",":")):
            full = read_json(archive/item["full"][-1])["result"]
            short = read_json(archive/item["short"][-1])["result"] if item["short"] else None
            curves.append((scoring_curve(full,short),name,color,ls))
        for data,name,color,ls in curves:
            outlet.plot(np.array(data["time_s"])/60,np.array(data["outlet_concentration_kg_m3"])/.15,
                        color=color,ls=ls,label=name,lw=1.7)
        outlet.set(xlim=(0,300),ylim=(0,1.02),xlabel="Time (min)")
        sink.grid(alpha=.16)
        outlet.grid(alpha=.16)
    axes[0,0].set_ylabel("Tracer captured in a period (%)")
    axes[1,0].set_ylabel("Outlet / inlet phenol concentration")
    axes[0,0].legend(loc="upper left",fontsize=7.4,frameon=False)
    handles,labels = axes[1,1].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=2,fontsize=8,frameon=False)
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"wall_uptake_figure7.png",dpi=180)
    plt.close(fig)
    write_json(HERE/"reference_summary.json",dict(summary=str(path.relative_to(ROOT)),sha256=sha(path),
        figure7=record["figure7"],uptake={size:{k:r[k] for k in ("uptake","step_changes","step_passed","coefficients")}
                                       for size,r in model["sizes"].items()},
        numerical={size:r["review"] for size,r in frozen["cases"].items()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve()
    render(path,validate_record(path))


if __name__=="__main__":
    main()
