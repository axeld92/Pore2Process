"""Plot the N128-to-N160 wall tracer comparison."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .native import HERE


def render(path):
    record = read_json(path)
    result = record["analysis"]
    sections = np.arange(1,9)
    fig, axes = plt.subplots(2,2,figsize=(10,6.1),sharex=True,sharey="row",layout="constrained")
    metrics = (("survival","All wall-free","#146c69"),
               ("fast_survival","Initially fast wall-free","#594a83"),
               ("conditional_capture","Capture among survivors","#a35d32"))
    for col,(size,diameter) in enumerate((("d800",.80),("d450",.45))):
        top,bottom = axes[:,col]
        primary,half = size+"_primary",size+"_half_step"
        for grid,color in (("N128","#4b5e86"),("N160","#237675")):
            data = result[grid]
            top.plot(sections,100*np.array(data[primary]["observed"]["fast_survival"]),
                     marker="o",ms=4,lw=1.8,color=color,label=grid+" · 80 µs")
            top.plot(sections,100*np.array(data[half]["observed"]["fast_survival"]),
                     ls="none",marker="x",ms=5,color=color,label=grid+" · 40 µs")
        bottom.axhspan(-2,2,color="#a8c6b9",alpha=.22)
        bottom.axvspan(.7,2.5,color="#cfcfcf",alpha=.25)
        bottom.axhline(0,color="#888888",lw=.75)
        for key,label,color in metrics:
            for variant,style in ((primary,"-"),(half,"--")):
                old,new = (np.asarray(result[g][variant]["observed"][key]) for g in ("N128","N160"))
                bottom.plot(sections,100*(new-old),ls=style,marker="o",ms=3,lw=1.5,
                            color=color,label=label if variant==primary else None)
        top.set_title(f"{diameter:.2f} mm particle-matched case")
        bottom.set(xlabel="Number of porous sections traversed",xticks=sections,xlim=(.75,8.25))
        for ax in (top,bottom):
            ax.grid(alpha=.17)
    axes[0,0].set_ylabel("Initially fast tracer wall-free (%)")
    axes[1,0].set_ylabel("N160 − N128 (percentage points)")
    axes[0,0].legend(ncol=2,fontsize=8,frameon=False,loc="lower left")
    handles,labels = axes[1,0].get_legend_handles_labels()
    for style,label in (("-","80 µs"),("--","40 µs")):
        handles.append(plt.Line2D([],[],ls=style,color="#333333"))
        labels.append(label)
    fig.legend(handles,labels,ncol=5,fontsize=8,frameon=False,loc="outside lower center")
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"disordered_n128_n160_tracer.png",dpi=175)
    plt.close(fig)
    write_json(HERE/"reference_summary.json",dict(path=str(path.relative_to(ROOT)),
        sha256=sha(path),parent_summary_sha256=record["parent_summary_sha256"],analysis=result))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path)
    render(parser.parse_args().summary.resolve())
