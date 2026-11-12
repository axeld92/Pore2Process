"""Plot fixed-inventory midpoint comparisons for both particle sizes."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

from .markers import HERE, OUTPUT


def render(path, record):
    fig, axes = plt.subplots(1,2,figsize=(10,3.55),layout="constrained",sharey=True)
    laws = (("none","No exchange","o","#8f572c"),
            ("synthetic","Original synthetic law","s","#237b78"),
            ("matched_Pe","Péclet-matched law","^","#634899"))
    for ax,(size,label) in zip(axes,(("d800","0.80 mm"),("d450","0.45 mm")),strict=True):
        data = record["analysis"][size]
        ax.axvline(data["observed_t50_s"]/60,color="#252525",linestyle="--",lw=1,label="Figure 7 midpoint")
        for law, name, marker, color in laws:
            values = [data["scenarios"][a][law]["t50_s"]/60 for a in ("flow","void","slow")]
            shift = {"none":-.17,"synthetic":0.,"matched_Pe":.17}[law]
            ax.scatter(values,[j+shift for j in range(3)],label=name,s=42,marker=marker,color=color)
        ax.set(title=label,xlabel="Time to 50% breakthrough (min)",
               yticks=[0,1,2],yticklabels=["Per unit flow","Per void volume","Slow-path enrichment"],
               ylim=(-.45,2.45))
        ax.grid(axis="x",alpha=.2)
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=4,fontsize=8,frameon=False)
    folder = HERE/"figures"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/"midpoint_attribution.png",dpi=180)
    plt.close(fig)
    base = path.parent
    compact = dict(summary=str(path.relative_to(ROOT)),sha256=sha(path),analysis=record["analysis"],
                   allocations=read_json(base/"allocation_summary.json"))
    write_json(HERE/"reference_summary.json",compact)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path)
    args = parser.parse_args()
    path = (args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"]).resolve()
    render(path,validate_record(path))


if __name__=="__main__":
    main()
