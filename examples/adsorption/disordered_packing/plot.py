"""Plot the same-grid ordered and disordered eight-cell tracer comparisons."""
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
    axes_keys = (("survival","All tracer wall-free (%)"),
                 ("fast_survival","Initially fast tracer wall-free (%)"))
    fig, axes = plt.subplots(2,2,figsize=(10,6.3),sharex=True,sharey="row",layout="constrained")
    x = np.arange(1,9)
    for j,(size,diameter) in enumerate((("d800",.80),("d450",.45))):
        for i,(key,ylabel) in enumerate(axes_keys):
            ax = axes[i,j]
            for name, label, color in (("ordered_N64","Ordered","#4b5e86"),
                                       ("disordered_N64","Disordered","#237675")):
                group = result["tracer"][name]
                ax.plot(x,100*np.array(group[size+"_primary"]["observed"][key]),
                        color=color,lw=2.,marker="o",ms=3.5,label=label+" · 80 µs")
                ax.plot(x,100*np.array(group[size+"_half_step"]["observed"][key]),
                        color=color,ls="none",marker="x",ms=5.,label=label+" · 40 µs")
            ax.set(xticks=x,xlim=(.75,8.25))
            ax.grid(alpha=.19)
            if j==0:
                ax.set_ylabel(ylabel)
        axes[0,j].set_title(f"{diameter:.2f} mm particle-matched Péclet number")
        axes[1,j].set_xlabel("Number of porous sections traversed")
    handles,labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=4,fontsize=9,frameon=False,loc="outside lower center")
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"ordered_disordered_tracer.png",dpi=175)
    plt.close(fig)
    write_json(HERE/"reference_summary.json",dict(path=str(path.relative_to(ROOT)),sha256=sha(path),
                                                carbon_mass_kg=record["carbon_mass_kg"],analysis=result))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path)
    render(p.parse_args().summary.resolve())
