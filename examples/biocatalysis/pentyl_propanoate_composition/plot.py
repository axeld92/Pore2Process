"""Plot the batch evidence, complete-condition prediction and bed diagnostics."""
from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import array,read_rows
from examples.biocatalysis.pentyl_propanoate_composition.run import HERE,PLAN,read_batch
from pore2process_demo.io import read_json


def plot(reference,output):
    summary=read_json(reference/"summary.json")
    pred_path=reference/"selected_batch_predictions.csv"
    if not pred_path.exists():
        pred_path=reference/"batch_predictions.csv"
    pred=read_rows(pred_path)
    groups=read_batch(read_json(PLAN))
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    blue,orange,red,gray,green="#17628C","#C26228","#A83A3A","#52687A","#407A50"
    ax=axs[0,0]
    for assay,color,mark,label in ((14,blue,"o","1:1, 60 °C"),(15,orange,"s","1:1, 70 °C"),
            (16,red,"^","1:1, 80 °C"),(12,gray,"D","2:1, 80 °C")):
        data=groups[assay]
        ax.plot(array(data,"time_min"),100*array(data,"printed"),marker=mark,ms=4,color=color,lw=1.4,label=label)
    ax.set(title="A  |  Equimolar feed loses rate when heated",xlabel="Elapsed time (min)",ylabel="Acid conversion (%)",
           ylim=(0,100))
    ax.legend(loc="best",fontsize=9)
    ax=axs[0,1]
    for assay,color,feed in ((12,blue,"2:1"),(16,orange,"1:1")):
        data=groups[assay]
        ax.plot(array(data,"time_min"),100*array(data,"printed"),"o",ms=4,color=color,label=f"Observed {feed}")
        for law,style,lbl in (("acid_decay","--","Complete loss"),("protected_residual_decay","-","Retained activity + protection")):
            selected=[r for r in pred if r["variant"]=="primary" and r["case"]=="temperature_80C"
                      and r["law"]==law and int(r["assay"])==assay]
            ax.plot(np.r_[0,array(selected,"time_min")],100*np.r_[0,array(selected,"predicted_conversion")],
                    color=color,lw=1.7,ls=style,label=f"{feed}: {lbl}")
    ax.set(title="B  |  Both 80 °C experiments held out",xlabel="Elapsed time (min)",ylabel="Acid conversion (%)",
           ylim=(0,100))
    ax.legend(fontsize=7.2,ncol=2,loc="best")
    ax=axs[1,0]
    laws=("ratio_reference","acid","acid_decay","protected_decay","residual_decay","protected_residual_decay")
    labels=("Ratio\nreference","Acid\ninhibition","Complete\nloss","Protected\nloss","Retained\nactivity","Protected\nretained")
    observed=[summary["batch"]["all_data"]["fit"][law]["equal_assay_rmse_pp"] for law in laws]
    withheld=[summary["batch"]["paired_scores"][law]["overall_equal_assay_rmse_pp"] for law in laws]
    idx=np.arange(len(laws))
    ax.bar(idx-.19,observed,.38,color=gray,label="Fit to all batch curves")
    ax.bar(idx+.19,withheld,.38,color=orange,label="Paired whole-curve holdout")
    ax.set_xticks(idx,labels,fontsize=8)
    ax.set(title="C  |  Fitting versus prediction",ylabel="Acid conversion RMSE (pp)",ylim=(0,1.18*max(withheld)))
    ax.legend(fontsize=9)
    ax=axs[1,1]
    cases=("prior_central_film_reference","all_batch_reversible_shared","active_batch_reversible_shared",
           "all_batch_protected_residual_decay_fully_active","all_batch_protected_residual_decay_long_time_floor")
    labels=("Prior\nfilm","Shared rate\n(all batch)","Shared rate\n(active batch)",
            "Protected\nfully active","Protected\nlong-time floor")
    rows={r["case"]:r for r in summary["sixty_C_three_ml_min"]}
    values=[100*rows[k]["predicted_conversion"] for k in cases]
    ax.bar(np.arange(len(cases)),values,color=[gray,orange,blue,green,green],width=.72)
    measured=100*rows[cases[0]]["observed_conversion"]
    ax.axhline(measured,ls=":",color=red,lw=2,label=f"Measured {measured:.1f}%")
    ax.set_xticks(np.arange(len(cases)),labels,fontsize=8)
    ax.set(title="D  |  Frozen 6 g bed scenarios — 60 °C, 3 mL/min",
           ylabel="Acid conversion (%)",ylim=(0,100))
    ax.legend(fontsize=9)
    for ax in axs.flat:
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",alpha=.16);ax.set_axisbelow(True)
    fig.suptitle("Separate local composition from catalyst exposure history",fontsize=16,fontweight="bold")
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=170,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--reference",type=Path,default=HERE/"reference")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    plot(args.reference,args.output)
