"""Plot the independent high-loading audit with physical axes."""
from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import array, read_rows
from examples.biocatalysis.pentyl_propanoate_loading_evidence.run import HERE, SOURCE
from pore2process_demo.io import read_json


def plot(reference, output):
    summary = read_json(reference/"summary.json")
    rows = read_rows(reference/"matched_progress.csv")
    fig, axes = plt.subplots(1,3,figsize=(13.5,4.35),layout="constrained")
    blue,orange,gray = "#17628C","#C26228","#52687A"
    ax=axes[0]
    palmitate = read_rows(SOURCE/"richetti2010/loading_progress.csv")
    for w,color,label in ((10,blue,"10% of substrate mass"),(20,orange,"20% of substrate mass")):
        x=array(palmitate,"time_min")*w/10
        ax.plot(x,array(palmitate,f"conversion_pct_at_{w}"),"o-",color=color,label=label)
    ax.axvline(30,ls=":",color=gray)
    ax.annotate("~46 point gap",xy=(30,72),xytext=(35,66),
                arrowprops={"arrowstyle":"-","color":gray},fontsize=10,color=gray)
    ax.set(xlim=(0,62),ylim=(0,104),xlabel="Time × (loading / 10%) (min)",
           ylabel="Acid conversion (%)",title="A  |  Palmitate: same catalyst exposure")
    ax.legend(loc="lower right",fontsize=8)
    ax=axes[1]
    for ratio,color,label in ((2,blue,"2:1, 25.4 → 50.9 g/L"),(1,orange,"1:1, 40.3 → 80.5 g/L")):
        selected=[r for r in rows if r["series"]==f"lorenzo2026_{ratio}to1" and float(r["loading_pct_reported"])==20]
        x=100*array(selected,"target_conversion")
        a=100*array(selected,"relative_activity_per_gram")
        lo=100*array(selected,"relative_activity_reading_low")
        hi=100*array(selected,"relative_activity_reading_high")
        ax.fill_between(x,lo,hi,color=color,alpha=.12)
        ax.plot(x,a,"o-",color=color,label=label)
    ax.axhline(100,ls=":",color=gray,label="Proportional")
    ax.set(xlabel="Matched acid conversion (%)",ylabel="High-loading activity per gram\n(relative to 10% acid-mass case, %)",
           ylim=(65,125),title="B  |  Pentyl propanoate: rate changes")
    ax.legend(loc="lower right",fontsize=8)
    ax=axes[2]
    keys=("film_mass_proportional","film_joint_pentyl_power_void",
          "film_external_only_power_mass_ratio_proxy","film_external_only_crowding_mass_ratio_proxy")
    names=["Proportional\n(reference)","Pentyl joint\npower","Other esters\npower","Other esters\ncrowding"]
    values=[summary["bed_comparison"][key]["rmse_percentage_points"] for key in keys]
    bars=ax.bar(np.arange(4),values,color=[gray,blue,orange,orange],width=.68)
    ax.set_xticks(np.arange(4),names,fontsize=8)
    ax.bar_label(bars,fmt="%.1f",padding=3,fontsize=9)
    ax.set(ylabel="Reactor prediction RMSE (percentage points)",ylim=(0,49),
           title="C  |  Transfer to fixed 6 g bed")
    for ax in axes:
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.16)
        ax.set_axisbelow(True)
    fig.suptitle("High-loading limitations exist, but do not transfer as one rate law",fontsize=15,fontweight="semibold")
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=175,bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--reference",type=Path,default=HERE/"reference")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    plot(args.reference,args.output)
