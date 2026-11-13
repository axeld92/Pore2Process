"""Plot the signed terms in the held-out variance audit."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import read_json

from .variance_audit import HERE


def main():
    record = read_json(HERE/"variance_reference.json")
    fig, axes = plt.subplots(2,2,figsize=(9.4,5.8),sharex=True,sharey=True,layout="constrained")
    terms = (("exchange","Path exchange","#237675"),("local_delay","Local transit times","#9b5433"),
             ("memory","Travel-time memory","#6d5088"),("propagation","Prior error carried forward","#757575"))
    for row,(size,diameter) in enumerate((("d800",.80),("d450",.45))):
        for col,(variant,step) in enumerate((("primary",80),("half_step",40))):
            ax = axes[row,col]
            rows = record["cases"][f"{size}_{variant}"]
            section = [item["section"] for item in rows]
            scale = 100/np.array([item["measured_variance_s2"] for item in rows])
            for key,label,color in terms:
                values = [item["cumulative_s2"][key] for item in rows]
                ax.plot(section,scale*values,color=color,label=label,lw=1.7)
            ax.plot(section,scale*[item["signed_error_s2"] for item in rows],color="#202633",
                    label="Observed − frozen variance",ls="--",marker="o",ms=3,lw=1.7)
            ax.axhline(0.,color="#aaa",lw=.8)
            ax.grid(alpha=.16)
            ax.set_title(f"{diameter:.2f} mm · {step} µs")
            ax.set(xticks=section,xlim=(2.85,8.15))
    fig.supylabel("Cumulative difference (% of observed variance)")
    for ax in axes[-1,:]:
        ax.set_xlabel("Held-out section")
    handles,labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=3,fontsize=8.5,frameon=False,loc="outside lower center")
    fig.savefig(HERE/"figures/variance_attribution.png",dpi=175)
    plt.close(fig)


if __name__=="__main__":
    main()
