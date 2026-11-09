"""Scientific comparison plots; uncertainty is unavailable, so no invented bars."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_results(result,output):
    styles={"calibration":("o","#16769a","Fitted: AKP-30 temperature series"),
            "primary_holdout":("s","#cf6712","Held out: AKP-30 hold-time series"),
            "transfer":("^","#6d5399","Transfer: AKP-15"),
            "out_of_domain":("X","#777777","Outside fit domain: CR-6 / AKP-50")}
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.9),layout="constrained")
    a,b=axes
    values=[c[key]/1e-16 for c in result["cases"] for key in ("measured_permeability_m2","predicted_permeability_m2")]
    lo,hi=min(values)*.65,max(values)*1.4
    xx=np.geomspace(lo,hi,100)
    a.fill_between(xx,.75*xx,1.25*xx,color="#dddddd",alpha=.6,label="±25% project screen",zorder=0)
    offsets={3:(-24,10),6:(18,-1),7:(-3,-18)}
    for role,(marker,color,label) in styles.items():
        group=[c for c in result["cases"] if c["role"]==role]
        x=np.array([c["measured_permeability_m2"] for c in group])/1e-16
        y=np.array([c["predicted_permeability_m2"] for c in group])/1e-16
        a.scatter(x,y,marker=marker,color=color,label=label,s=65)
        for c,xx,yy in zip(group,x,y,strict=True):
            a.annotate(str(c["table_row"]),(xx,yy),xytext=offsets.get(c["table_row"],(5,5)),textcoords="offset points",fontsize=9,
                       arrowprops={"arrowstyle":"-","color":"#777777","lw":.6} if c["table_row"] in offsets else None)
    a.plot([lo,hi],[lo,hi],color="black",lw=1,label="Prediction = measurement")
    a.set(xscale="log",yscale="log",xlim=(lo,hi),ylim=(lo,hi),xlabel="Measured support k (10⁻¹⁶ m²)",ylabel="Predicted support k (10⁻¹⁶ m²)",title="One fitted multiplier; labels are Table 2 rows")
    a.legend(fontsize=7,loc="upper left")
    primary=[c for c in result["cases"] if c["role"]=="primary_holdout"]
    positions=np.arange(len(primary))
    b.axhspan(-25,25,color="#dddddd",alpha=.6)
    b.bar(positions-.17,[100*c["baseline_relative_error"] for c in primary],width=.34,color="#aaaaaa",label="Unfitted Carman–Kozeny")
    b.bar(positions+.17,[100*c["relative_error"] for c in primary],width=.34,color="#cf6712",label="Calibrated closure")
    b.axhline(0,color="black",lw=1)
    b.axhline(25,color="#777777",ls="--",lw=1)
    b.axhline(-25,color="#777777",ls="--",lw=1)
    b.set(xticks=positions,xticklabels=[f"Row {c['table_row']}" for c in primary],ylabel="(Predicted / measured − 1) (%)",title="Held-out predictions; shaded band is project screen")
    b.legend(fontsize=8)
    fig.suptitle("Biesheuvel & Verweij support measurements — retrospective same-study check",fontsize=12)
    fig.savefig(output/"permeability_validation.png",dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.8,4.8),layout="constrained")
    x=np.array([c["selective_to_measured_baseline_support_ratio"] for c in result["conditional_process"]])
    for key,label,style in (("measured_support_k","Measured support k; assumed fixed layer","o-"),
                            ("calibrated_support_k","Fitted support k; same assumed fixed layer","s--")):
        y=[100*c[key]["fractional_shaft_energy_saving"] for c in result["conditional_process"]]
        ax.plot(x,y,style,label=label)
    ax.set(xscale="log",xlabel="Fixed selective resistance / measured baseline support resistance (1)",
           ylabel="Conditional shaft-energy saving (%)",title="AKP-30 row 3 → row 2; 100 LMH, nominal 2 mm supports")
    ax.legend(fontsize=9)
    ax.text(.02,.02,"Illustrative process calculation; not measured whole-membrane energy",transform=ax.transAxes,fontsize=9)
    ax.grid(alpha=.2)
    fig.savefig(output/"conditional_energy.png",dpi=180)
    plt.close(fig)
