"""Plot held-out variance errors from the compact N160 reference."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import read_json
from .run import HERE


def main():
    record = read_json(HERE/"reference_summary.json")
    sections = np.arange(3,9)
    fig,axes = plt.subplots(2,2,figsize=(9,5.5),sharex=True,sharey=True,layout="constrained")
    for row,(size,diameter) in enumerate((("d800",.80),("d450",.45))):
        for col,(variant,step) in enumerate((("primary",80),("half_step",40))):
            ax = axes[row,col]
            case = record["analysis"]["cases"][f"{size}_{variant}"]
            gate = 100*record["analysis"]["limits"]["variance"]
            ax.axhspan(-gate,gate,color="#eef4ef",label="±10% limit")
            ax.axhline(0,color="#85928b",lw=.7)
            a,b = [100*np.array(case[k]["variance_relative_error"]) for k in ("probe_first","probe_second")]
            ax.fill_between(sections,a,b,color="#a3cfcb",alpha=.6,label="Two probe halves")
            for key,label,color,marker in (
                ("stationary","Repeated rule","#9b5433","s"),
                ("sectional","Section-specific rules","#146a65","o"),
            ):
                ax.plot(sections,100*np.array(case[key]["variance_relative_error"]),color=color,
                        label=label,marker=marker,ms=3.5,lw=1.6)
            ax.set(title=f"{diameter:.2f} mm · {step} µs",xticks=sections,xlim=(2.85,8.15))
            ax.grid(axis="x",alpha=.12)
    for ax in axes[-1,:]:
        ax.set_xlabel("Held-out section")
    fig.supylabel("Arrival-time variance error (%)")
    handles,labels = axes[0,0].get_legend_handles_labels()
    order = (2,3,1,0)
    fig.legend([handles[i] for i in order],[labels[i] for i in order],ncol=4,
               fontsize=8.5,frameon=False,loc="outside lower center")
    path = HERE/"figures/variance_holdout.png"
    path.parent.mkdir(exist_ok=True)
    fig.savefig(path,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    main()
