"""Render the N160 holdout and save the compact audited result."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .run import HERE


def render(path):
    record = validate_record(path)
    frozen = read_json(path.parent/"frozen_prediction.json")
    cases = record["analysis"]["cases"]
    cells = np.arange(3,9)
    fig, axes = plt.subplots(2,2,figsize=(9.5,5.9),sharex=True,sharey="row",layout="constrained")
    for j,(size,diameter) in enumerate((("d800",.80),("d450",.45))):
        top,bottom = axes[:,j]
        for variant,style in (("primary","-"),("half_step","--")):
            case = cases[f"{size}_{variant}"]
            pred,obs = case["pair"]["profiles"],case["observed"]
            for key,label,color in (("survival","All wall-free","#237675"),
                                   ("fast_survival","Initially fast wall-free","#594a83"),
                                   ("conditional_capture","Capture among survivors","#a35d32")):
                diff = np.asarray(pred[key][2:])-obs[key][2:]
                top.plot(cells,100*diff,color=color,ls=style,marker="o",ms=3,
                         label=label if variant=="primary" else None)
            for key,label,color in (("mean_arrival_s","Mean arrival time","#4b5e86"),
                                   ("variance_s2","Arrival-time variance","#ad4f3f")):
                ratio = np.asarray(pred[key][2:])/obs[key][2:]-1
                bottom.plot(cells,100*ratio,color=color,ls=style,marker="o",ms=3,
                            label=label if variant=="primary" else None)
        top.axhspan(-2,2,color="#a8c6b9",alpha=.23)
        bottom.axhspan(-10,10,color="#a8c6b9",alpha=.23)
        for y in (-5,5):
            bottom.axhline(y,color="#4b5e86",alpha=.45,lw=.8)
        for ax in (top,bottom):
            ax.axhline(0,color="#777777",lw=.7)
            ax.grid(alpha=.15)
        top.set_title(f"{diameter:.2f} mm particle-matched case")
        bottom.set(xlabel="Held-out section",xticks=cells,xlim=(2.8,8.2))
    axes[0,0].set_ylabel("Wall-contact error (percentage points)")
    axes[1,0].set_ylabel("Inert timing error (%)")
    axes[0,0].legend(fontsize=8,frameon=False,loc="upper left")
    handles,labels = axes[1,0].get_legend_handles_labels()
    for style,label in (("-","80 µs"),("--","40 µs")):
        handles.append(plt.Line2D([],[],ls=style,color="#333333"))
        labels.append(label)
    fig.legend(handles,labels,ncol=4,fontsize=8,frameon=False,loc="outside lower center")
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"two_section_holdout.png",dpi=175)
    plt.close(fig)
    coefficients = {key:{k:value[k] for k in ("Pe","dt_s","training_sha256","counts","clear",
                    "delay","delay2","states","flux","fast_share")} for key,value in frozen["cases"].items()}
    write_json(HERE/"reference_summary.json",dict(path=str(path.relative_to(ROOT)),sha256=sha(path),
        frozen_sha256=record["frozen_prediction_sha256"],parent_sha256=record["parent_sha256"],
        training_coefficients=coefficients,analysis=record["analysis"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path)
    render(parser.parse_args().summary.resolve())
