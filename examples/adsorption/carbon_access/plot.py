"""Plot the midpoint and retained mass for the paired carbon-access scenarios."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import read_json
from .run import HERE


def plot():
    record = read_json(HERE/"reference_summary.json")
    fig,axes = plt.subplots(2,2,figsize=(9.0,5.5),sharex=True,sharey="row",layout="constrained")
    positions = np.arange(3)
    for col,(size,diameter) in enumerate((("d800",.80),("d450",.45))):
        case = record["figure7"][size]
        full = case["cases"]["full"]
        top,bottom = axes[:,col]
        m0 = case["figure7_midpoint_s"]/60
        r0 = case["figure7_retained_kg"]*1e6
        top.axhspan(.8*m0,1.2*m0,color="#e5e8eb",label="Figure 7 comparison range")
        bottom.axhspan(*np.array(case["reading_scenarios_kg"])*1e6,color="#e5e8eb")
        for ax,value,control in ((top,m0,full["midpoint_s"]/60),(bottom,r0,full["snapshot"]["retained_kg"]*1e6)):
            ax.axhline(control,color="#92989d",ls="--",label="Full direct access")
            ax.axhline(value,color="#34405a",ls=":",label="Figure 7 value")
            ax.grid(axis="y",alpha=.15)
            ax.set_xlim(-.25,2.25)
        for frac,label,color,symbol in ((50,"50% direct access","#136d67","o"),
                                         (75,"75% direct access","#975134","s")):
            rows = [case["cases"][f"f{frac}_{name}"] for name in ("isolated","slow01","slow10")]
            top.plot(positions,[r["midpoint_s"]/60 for r in rows],color=color,marker=symbol,label=label)
            bottom.plot(positions,[1e6*r["snapshot"]["retained_kg"] for r in rows],color=color,marker=symbol)
        top.set_title(f"{diameter:.2f} mm")
        bottom.set(xticks=positions,xticklabels=("Isolated","1%","10%"),
                   xlabel="Remaining carbon’s conductance relative to direct film")
    axes[0,0].set_ylabel("Time to 50% of feed (min)")
    axes[1,0].set_ylabel("Retained phenol (mg)")
    handles,labels = axes[0,0].get_legend_handles_labels()
    ii = (3,4,1,2,0)
    fig.legend([handles[i] for i in ii],[labels[i] for i in ii],loc="outside lower center",ncol=3,
               frameon=False,fontsize=9)
    dest = HERE/"figures/carbon_access.png"
    dest.parent.mkdir(exist_ok=True)
    fig.savefig(dest,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    plot()
