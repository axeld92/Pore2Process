"""Plot the N96-to-N128 wall-free tracer and the two-point grid screen."""
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
        for grid,color in (("N96","#4b5e86"),("N128","#237675")):
            data = result[grid]
            top.plot(sections,100*np.array(data[primary]["observed"]["fast_survival"]),
                     marker="o",ms=4,lw=1.8,color=color,label=grid+" · 80 µs")
            top.plot(sections,100*np.array(data[half]["observed"]["fast_survival"]),
                     ls="none",marker="x",ms=5,color=color,label=grid+" · 40 µs")
        bottom.axhspan(-2,2,color="#a8c6b9",alpha=.22)
        bottom.axvspan(.7,2.5,color="#cfcfcf",alpha=.25)
        bottom.axhline(0,color="#888888",lw=.75)
        for key,label,color in metrics:
            old,new = (result[grid][primary]["observed"][key] for grid in ("N96","N128"))
            bottom.plot(sections,100*(np.array(new)-old),marker="o",ms=3.4,
                        lw=1.7,color=color,label=label)
        top.set_title(f"{diameter:.2f} mm particle-matched case")
        bottom.set(xlabel="Number of porous sections traversed",xticks=sections,xlim=(.75,8.25))
        for ax in (top,bottom):
            ax.grid(alpha=.17)
    axes[0,0].set_ylabel("Initially fast tracer wall-free (%)")
    axes[1,0].set_ylabel("N128 − N96 at 80 µs (percentage points)")
    axes[0,0].legend(ncol=2,fontsize=8,frameon=False,loc="lower left")
    handles,labels = axes[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=3,fontsize=9,frameon=False,loc="outside lower center")
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"disordered_n96_n128_tracer.png",dpi=175)
    plt.close(fig)
    write_json(HERE/"reference_summary.json",dict(path=str(path.relative_to(ROOT)),
        sha256=sha(path),parent_summary_sha256=record["parent_summary_sha256"],analysis=result))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path)
    render(p.parse_args().summary.resolve())
