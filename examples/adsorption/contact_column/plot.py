"""Plot the two unfitted contact-based column predictions and source markers."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.alhamed.study import observations
from examples.adsorption.exchange_column.run import CASES
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .inputs import ATTRIBUTION
from .run import HERE, OUTPUT, scoring_curve

STYLES = {
    "advection_rank": ("Contact allocation, advection ranks","#237b78","-"),
    "realized_rank": ("Contact allocation, diffusive ranks","#694499","-."),
}


def render(path,record,frozen,model):
    base = path.parent
    fig, axes = plt.subplots(2,2,figsize=(10,6.2),layout="constrained",sharey="row",
                             height_ratios=(1.45,1))
    prior = read_json(ATTRIBUTION/"frozen_predictions.json")["cases"]
    for column,(size,(series,diameter)) in enumerate(CASES.items()):
        a,b = axes[:,column]
        t,y = observations(series)
        a.scatter(t/60,y,facecolors="none",edgecolors="#252525",s=19,label="Figure 7 markers",zorder=4)
        mask = t<=60
        b.scatter(t[mask],y[mask],facecolors="none",edgecolors="#252525",s=19,zorder=4)
        previous = read_json(ATTRIBUTION/prior[f"{size}_void_matched_Pe"]["source"])
        cases = [(previous["result"],"Void-proportional carbon","#95623a","--")]
        for rank,item in frozen["cases"][size].items():
            full = read_json(base/item["selected"])["result"]
            early = read_json(base/item["startup_files"][-1])["result"] if item["startup_files"] else None
            cases.append((scoring_curve(full,early),*STYLES[rank]))
        for result,label,color,style in cases:
            time = np.array(result["time_s"])
            value = np.array(result["outlet_concentration_kg_m3"])/previous["inputs"]["inlet_concentration_kg_m3"]
            a.plot(time/60,value,color=color,ls=style,lw=1.65,label=label)
            b.plot(time,value,color=color,ls=style,lw=1.65)
        a.set(xlim=(0,300),ylim=(0,1.02),title=f"{1000*diameter:g} mm particles",xlabel="Time (min)")
        b.set(xlim=(0,60),ylim=(0,1.02),xlabel="First minute (s)")
        for ax in (a,b):
            ax.grid(alpha=.17)
    axes[0,0].set_ylabel("Outlet / inlet phenol concentration")
    axes[1,0].set_ylabel("Outlet / inlet phenol concentration")
    handles, labels = axes[0,1].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=2,fontsize=8,frameon=False)
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"figure7_contact.png",dpi=180)
    plt.close(fig)
    cases = {}
    for size,settings in model["sizes"].items():
        cases[size] = dict(allocations=settings["allocation_summary"],figure7=record["figure7"][size],
                           numerical_checks={rank:{k:v for k,v in case.items() if k!="files"}
                                             for rank,case in frozen["cases"][size].items()})
    write_json(HERE/"reference_summary.json",dict(summary=str(path.relative_to(ROOT)),sha256=sha(path),
               dry_carbon_mass_kg=model["dry_carbon_mass_kg"],cases=cases))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve()
    render(path,validate_record(path),read_json(path.parent/"frozen_predictions.json"),
           read_json(path.parent/"model.json"))


if __name__=="__main__":
    main()
