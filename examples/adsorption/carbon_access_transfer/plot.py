"""Compare retained mass against inlet dose on the reported dry-mass basis."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import read_json, sha
from . import source
from .run import HERE, OUTPUT, ROOT


def plot(path):
    record = read_json(path)
    base = path.parent
    if sha(base/"frozen.json")!=record["frozen_sha256"]:
        raise ValueError("Frozen curve manifest changed")
    setup,frozen = [read_json(base/name) for name in ("plan.json","frozen.json")]
    markers = source.markers()
    variants = (
        ("full","Full access","#777c82","--"),
        ("f50_isolated","Half isolated","#9b5431",":"),
        ("f50_slow01","Half at 1% reference conductance","#5c5285","-."),
        ("f50_slow10","Half at 10% reference conductance","#176e67","-"),
    )
    dest = HERE/"figures"
    dest.mkdir(exist_ok=True)
    for figure in (5,4):
        shape = (2,2) if figure==5 else (1,3)
        common = 1e3*record["score"]["common_dose_kg_kg"][str(figure)]
        fig,axes = plt.subplots(*shape,figsize=(9.8,6.2 if figure==5 else 4.1),
                                sharex=True,sharey=True,layout="constrained")
        items = [(key,m) for key,m in markers.items() if m["figure"]==figure]
        for ax,(key,m) in zip(axes.flat,items,strict=True):
            case = m["case"]
            cond = setup["cases"][case]
            p = cond["parameters"]
            scale = 1e3/p["dry_mass_kg"]
            feed = p["flow_m3_s"]*p["inlet_concentration_kg_m3"]
            stop = m["time"][-1]
            for name,label,color,style in variants:
                file = frozen["cases"][case][name]["selected"]
                if sha(base/file)!=frozen["curves_sha256"][file]:
                    raise ValueError("Frozen column curve changed")
                r = read_json(base/file)["result"]
                t = np.array(r["time_s"])
                t_view = np.r_[t[t<stop],stop]
                retained = np.array(r["adsorbed_inventory_kg"])+r["liquid_inventory_kg"]
                ax.plot(t_view*feed*scale,np.interp(t_view,t,retained)*scale,color=color,ls=style,lw=1.6,label=label)
            t,frac = np.r_[0.,m["time"]],np.r_[0.,m["fraction"]]
            retained = np.cumsum(np.diff(t)*(2-frac[1:]-frac[:-1])*feed*scale/2)
            ax.plot(t[1:]*feed*scale,retained,"o",ms=3,mec="#292d34",mfc="white",label="Published markers")
            ax.axvline(common,color="#929ca8",lw=.8,label="Common feed dose")
            label = f'{cond["bed_cm"]} cm' if figure==5 else f'{cond["flow_mL_min"]} mL/min'
            ax.set_title(label)
            ax.grid(alpha=.15)
        for ax in np.asarray(axes).reshape(shape)[-1]:
            ax.set_xlabel("Phenol fed (mg/g carbon)")
        for ax in np.asarray(axes).reshape(shape)[:,0]:
            ax.set_ylabel("Phenol retained (mg/g carbon)")
        handles,labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles,labels,ncol=2,fontsize=8.7,frameon=False,loc="outside lower center")
        fig.savefig(dest/f"figure{figure}_retention.png",dpi=170)
        plt.close(fig)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path,nargs="?")
    args = p.parse_args()
    plot(args.summary or ROOT/read_json(OUTPUT/"latest.json")["summary"])
