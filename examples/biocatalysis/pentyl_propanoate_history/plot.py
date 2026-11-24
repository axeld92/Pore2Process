"""Chemical-engineering figures: conversion over time and stationary catalyst state."""
from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import array,read_rows
from examples.biocatalysis.pentyl_propanoate_history.run import HERE

BLUE,ORANGE,TEAL,GRAY="#17628C","#C26228","#337A65","#52687A"
LAW_STYLES=( ("acid_decay","Unprotected damage",BLUE,"--"),
             ("acid_both","Damage + recovery",ORANGE,"-"),
             ("protected_residual_decay","Alcohol protection + residual activity",TEAL,"-") )


def condition(rows,feed,T,Q):
    return [r for r in rows if int(r["feed_ratio"])==feed and int(r["temperature_C"])==T and float(r["flow_ml_min"])==Q]


def style(ax):
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y",alpha=.18);ax.set_axisbelow(True)


def plot(reference,out):
    history=[r for r in read_rows(reference/"history.csv") if r["scenario"]=="fresh_bulk"]
    comp=[r for r in read_rows(reference/"bed_comparison.csv") if r["scenario"]=="fresh_bulk"]
    profile=read_rows(reference/"selected_axial_profiles.csv")
    out.mkdir(parents=True,exist_ok=True)
    fig,axs=plt.subplots(2,2,figsize=(12.4,8.2))
    for j,(feed,T,Q) in enumerate(((1,80,.7),(2,60,3.))):
        ax=axs[0,j]
        for law,label,color,ls in LAW_STYLES:
            r=sorted([r for r in condition(history,feed,T,Q) if r["law"]==law],key=lambda r:float(r["time_min"]))
            ax.plot(array(r,"time_min"),100*array(r,"predicted_conversion"),color=color,lw=2,ls=ls,marker="o",ms=3,label=label)
        row=condition(comp,feed,T,Q)[0]
        ax.axhline(100*float(row["observed_conversion"]),color="black",ls=":",lw=1.7,label="Reported conversion (age unknown)")
        ax.axvspan(30,180,color=GRAY,alpha=.06)
        ax.axvline(30,color=GRAY,ls=":",lw=1)
        ax.set(xlim=(0,180),xticks=[0,30,60,90,120,150,180],ylim=(0,100),xlabel="Time on stream (min)",
               title=f"{feed}:1 • {T} °C • {Q:g} mL/min",ylabel="Outlet acid conversion (%)" if j==0 else "")
        style(ax)
        ax=axs[1,j]
        law="acid_both" if j==0 else "protected_residual_decay"
        for age,color,ls in ((10,BLUE,"--"),(30,ORANGE,"-"),(180,TEAL,"-")):
            r=sorted([r for r in condition(profile,feed,T,Q) if r["law"]==law and float(r["time_min"])==age],
                     key=lambda r:float(r["axial_fraction"]))
            ax.plot(140*array(r,"axial_fraction"),100*array(r,"activity"),color=color,lw=2,ls=ls,label=f"{age} min")
        label="Damage + recovery" if j==0 else "Alcohol protection + residual activity"
        ax.set(xlim=(0,140),ylim=(0,105),xlabel="Axial distance from inlet (mm)",
               ylabel="Relative catalyst activity (%)" if j==0 else "",title=label)
        ax.legend(loc="lower right",fontsize=9);style(ax)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,.93),ncol=2,frameon=False,fontsize=10)
    fig.suptitle("Stationary catalyst history | Same 6 g bed; no reactor fitting",fontsize=17,fontweight="bold",y=.99)
    fig.text(.5,.015,"Shaded time window: 30–180 min. Horizontal markers are single reported conversions, not measured time traces.",ha="center",fontsize=9)
    fig.subplots_adjust(left=.075,right=.985,top=.80,bottom=.085,wspace=.18,hspace=.34)
    fig.savefig(out/"local_history.png",dpi=170,bbox_inches="tight");plt.close(fig)

    fig,axs=plt.subplots(2,4,figsize=(14.4,7.2),sharex=True,sharey=True)
    for i,feed in enumerate((1,2)):
        for j,Q in enumerate((.7,1.,2.,3.)):
            ax=axs[i,j]
            rows=sorted([r for r in comp if int(r["feed_ratio"])==feed and float(r["flow_ml_min"])==Q and
                         r["law"]=="protected_residual_decay"],key=lambda r:float(r["temperature_C"]))
            T=array(rows,"temperature_C")
            ax.fill_between(T,100*array(rows,"minimum_30to180_conversion"),100*array(rows,"maximum_30to180_conversion"),
                alpha=.16,color=TEAL,label="Protected model, 30–180 min span")
            ax.plot(T,100*array(rows,"previous_fully_active_steady_conversion"),ls=":",color=GRAY,lw=2,label="Protected model, fully active limit")
            for law,label,color,ls in LAW_STYLES:
                r=sorted([r for r in comp if int(r["feed_ratio"])==feed and float(r["flow_ml_min"])==Q and r["law"]==law],
                         key=lambda r:float(r["temperature_C"]))
                ax.plot(T,100*array(r,"predicted_30to180_sample_mean"),color=color,lw=1.9,ls=ls,label=label)
            ax.plot(T,100*array(rows,"observed_conversion"),"ko",ms=4.7,label="Held-out observations")
            ax.set(xlim=(37,83),ylim=(0,100),xticks=[40,50,60,70,80],title=f"{feed}:1  |  {Q:g} mL/min")
            if j==0:ax.set_ylabel("Acid conversion (%)")
            if i==1:ax.set_xlabel("Temperature (°C)")
            style(ax)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,.93),ncol=3,fontsize=9,frameon=False)
    fig.suptitle("All forty reactor points | Mean of six prespecified 30–180 min predictions",fontsize=16,fontweight="bold",y=.987)
    fig.subplots_adjust(left=.065,right=.99,top=.79,bottom=.105,wspace=.13,hspace=.24)
    fig.savefig(out/"dual_feed_history.png",dpi=160,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,default=HERE/"reference")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();plot(args.reference,args.output)
