"""Plot the full dual-feed experiment matrix and the two diagnostic contrasts."""
from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import array,read_rows
from examples.biocatalysis.pentyl_propanoate_equimolar.run import HERE

BLUE, TEAL, ORANGE, GRAY = "#17628C", "#337A65", "#C26228", "#52687A"
FLOWS = (.7,1.,2.,3.)


def line(rows,feed,Q,key,value):
    r=sorted([d for d in rows if int(d["feed_ratio"])==feed and float(d["flow_ml_min"])==Q
              and d[key]==value],key=lambda d:float(d["temperature_C"]))
    if len(r)!=5:
        raise ValueError(f"Expected five temperatures: {feed}, {Q}, {value}")
    return r


def style(ax):
    ax.set(xlim=(37,83),ylim=(0,100),xticks=[40,50,60,70,80])
    ax.spines[["top","right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y",alpha=.2)


def plot(source,sign_source,output):
    rows=read_rows(source/"bed_comparison.csv")
    clock=read_rows(source/"clock_comparison.csv")
    sign=read_rows(sign_source/"sign_bed_comparison.csv")
    output.mkdir(parents=True,exist_ok=True)
    # Full condition matrix, with the same axes for the two feeds.
    fig,axs=plt.subplots(2,4,figsize=(14.4,7.3),sharex=True,sharey=True)
    for i,feed in enumerate((1,2)):
        for j,Q in enumerate(FLOWS):
            ax=axs[i,j]
            for case,label,color,ls in (
                ("active_batch_reversible_shared","Shared active-domain fit",BLUE,"--"),
                ("all_batch_ratio_reference","Separate feed-rate fit",ORANGE,":"),
                ("all_batch_composition","Composition fit to all batches",TEAL,"-")):
                r=line(rows,feed,Q,"case",case)
                ax.plot(array(r,"temperature_C"),100*array(r,"predicted_conversion"),ls=ls,color=color,lw=2,label=label)
            ax.plot(array(r,"temperature_C"),100*array(r,"observed_conversion"),"ko",ms=4.5,label="Held-out bed measurements")
            ax.set_title(f"{feed}:1 feed  |  {Q:g} mL/min",fontsize=11,loc="left")
            if j==0: ax.set_ylabel("Acid conversion (%)")
            if i==1: ax.set_xlabel("Temperature (°C)")
            style(ax)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,.93),ncol=4,frameon=False,fontsize=10)
    fig.suptitle("Independent bed test • fixed 6 g catalyst • no reactor fitting",fontsize=17,fontweight="bold",y=.985)
    fig.subplots_adjust(left=.07,right=.99,top=.83,bottom=.10,wspace=.13,hspace=.27)
    fig.savefig(output/"dual_feed.png",dpi=160,bbox_inches="tight");plt.close(fig)

    # Two mechanistically distinct questions at the same flow rate.
    fig,axs=plt.subplots(2,2,figsize=(10.8,8.1),sharex=True,sharey=True)
    for j,feed in enumerate((1,2)):
        ax=axs[0,j]
        a=line(sign,feed,3.,"case","literal_printed")
        b=line(sign,feed,3.,"case","opposite_exponent_sign")
        T=array(a,"temperature_C")
        ax.plot(T,100*array(a,"author_predicted_conversion"),color=GRAY,ls=":",lw=2,label="Paper's calculated curve")
        ax.plot(T,100*array(a,"predicted_conversion"),color=BLUE,ls="--",lw=2,label="Literal printed sign, 6 g")
        ax.plot(T,100*array(b,"predicted_conversion"),color=ORANGE,lw=2,label="Opposite sign, 6 g")
        ax.plot(T,100*array(a,"observed_conversion"),"ko",ms=5,label="Measured")
        ax.set_title(f"{feed}:1  |  Published rate sign audit",loc="left")
        style(ax)
        ax=axs[1,j]
        full=line(rows,feed,3.,"case","all_batch_protected_residual_decay_fully_active")
        for age,col,alpha in ((180.,ORANGE,.17),(30.,BLUE,.21)):
            c=line([r for r in clock if float(r["catalyst_age_min"])==age],feed,3.,"law","protected_residual_decay")
            key="conversion_uniform_inlet_activity"; key2="conversion_uniform_equilibrium_activity"
            ax.fill_between(T,100*array(c,key),100*array(c,key2),color=col,alpha=alpha,label=f"{age:g} min uniform activity range")
            ax.plot(T,100*array(c,key),color=col,lw=1.2)
        ax.plot(T,100*array(full,"predicted_conversion"),ls="--",lw=1.7,color=TEAL,label="Fully active limit")
        ax.plot(T,100*array(full,"observed_conversion"),"ko",ms=5,label="Measured")
        ax.set_title(f"{feed}:1  |  Stationary catalyst-age scenarios",loc="left")
        ax.set_xlabel("Temperature (°C)")
        style(ax)
    for i in (0,1):
        axs[i,0].set_ylabel("Acid conversion (%)")
        axs[i,1].legend(loc="best",fontsize=8,frameon=True)
    fig.suptitle("Composition and catalyst activity at 3 mL/min",fontsize=16,fontweight="bold",y=.985)
    fig.text(.5,.012,"Shaded ranges are conditional exposure scenarios, not measurement uncertainty or fitted bed activity.",ha="center",fontsize=9)
    fig.subplots_adjust(left=.09,right=.98,top=.91,bottom=.09,wspace=.14,hspace=.24)
    fig.savefig(output/"source_and_age.png",dpi=170,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE/"reference")
    p.add_argument("--sign-source",type=Path,default=HERE/"reference")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();plot(args.source,args.sign_source,args.output)
