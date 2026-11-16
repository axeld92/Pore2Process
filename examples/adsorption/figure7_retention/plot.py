"""Plot cumulative phenol retention and the change in removal rate with loading."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .run import HERE, figure_curve, inputs, integrate


def plot():
    _,basis,markers,models,_ = inputs()
    feed = basis["feed_rate_kg_s"]
    fig,axes = plt.subplots(2,2,figsize=(9.3,5.9),sharey="row",layout="constrained")
    for col,(size,data) in enumerate(markers.items()):
        top,bottom = axes[:,col]
        time,fraction = figure_curve(data)
        dense = np.unique(np.r_[time,np.linspace(0.,time[-1],600)])
        accumulated = lambda t: np.array([1e6*integrate(time,fraction,feed,s) for s in t])
        for key,label,color,style in (
                ("mixed","Fully mixed joins","#9b5334","--"),
                ("preserved","Path-preserving joins","#136d67","-")):
            m = models[size][key]
            t = m["time_s"]
            window = np.r_[t[t<time[-1]],time[-1]]
            retained = 1e6*np.interp(window,t,m["retained_kg"])
            top.plot(window/60,retained,color=color,linestyle=style,label=label)
            bottom.plot(retained,60e6*feed*(1-np.interp(window,t,m["outlet_fraction"])),
                        color=color,linestyle=style)
        top.plot(dense/60,accumulated(dense),color="#343944",lw=1.)
        top.scatter(time[1:]/60,accumulated(time[1:]),s=14,facecolor="white",edgecolor="#343944",
                    linewidth=.9,zorder=3,label="Figure 7 integration")
        bottom.plot(accumulated(dense),60e6*feed*(1-np.interp(dense,time,fraction)),color="#343944",lw=1.)
        bottom.scatter(accumulated(time[1:]),60e6*feed*(1-fraction[1:]),s=14,
                       facecolor="white",edgecolor="#343944",linewidth=.9,zorder=3)
        top.axhline(1e6*basis["equilibrium_total_inventory_kg"],color="#999",ls=":",label="Full equilibrium inventory")
        top.set(title=f'{1000*data["diameter_m"]:.2f} mm',xlabel="Time (min)",xlim=(0,300),ylim=(0,800))
        bottom.set(xlabel="Retained phenol (mg)",xlim=(0,800),ylim=(0,7))
        for axis in (top,bottom):
            axis.grid(alpha=.16)
    axes[0,0].set_ylabel("Retained phenol (mg)")
    axes[1,0].set_ylabel("Net removal (mg/min)")
    handles,labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=2,frameon=False,fontsize=9)
    dest = HERE/"figures/figure7_retained_phenol.png"
    dest.parent.mkdir(exist_ok=True)
    fig.savefig(dest,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    plot()
