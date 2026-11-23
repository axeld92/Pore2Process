"""Plot the four measured progress curves without kinetic fitting."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from examples.biocatalysis.pentyl_propanoate.run import BATCH, array, read_rows
from examples.biocatalysis.pentyl_propanoate_loading.run import PLAN, batch_rows, matched_conversion
from pore2process_demo.io import read_json


def plot(path):
    plan = read_json(PLAN)
    groups = batch_rows(read_rows(BATCH), plan)
    _, matched = matched_conversion(groups, plan)
    reference = plan["batch_source"]["reference_loading_pct"]
    levels = plan["matched_conversion_levels"]
    figure = Figure(figsize=(13.2, 4.5), layout="constrained", dpi=160)
    FigureCanvasAgg(figure)
    actual, scaled, rate = figure.subplots(1,3)
    colors = ("#2478b4", "#279265", "#6636a4", "#d17b18")
    markers = ("o","^","s","D")
    for (loading,rows),color,marker in zip(groups.items(),colors,markers,strict=True):
        time, conversion = array(rows,"time_min"),array(rows,"printed")
        options = dict(color=color,marker=marker,markersize=4,linewidth=1.5,label=f"{loading:g}%")
        actual.plot(time,100*conversion,**options)
        scaled.plot(time*loading/reference,100*conversion,**options)
        ratio = [matched["printed"][str(X)]["cumulative_mass_specific_rate_relative_to_10pct"][str(loading)]
                 for X in levels]
        rate.plot(np.asarray(levels)*100,ratio,**options)
    actual.set(title="A  Measured batch curves",xlabel="Elapsed time (min)",ylabel="Acid conversion (%)",
               xlim=(0,330),ylim=(0,90))
    actual.legend(title="Catalyst / initial acid mass",loc="lower right",fontsize=8,title_fontsize=8)
    # Lowest loading's last observation determines support common to all.
    # Use the same y-axis and avoid extending its data by extrapolation.
    common_end = min(rows[-1]["time_min"]*w/reference for w,rows in groups.items())
    scaled.set(title="B  Time × relative catalyst mass",xlabel="Equivalent time at 10% loading (min)",
               xlim=(0,common_end*1.1),ylim=(0,90))
    scaled.axvline(common_end,color="#828282",linewidth=1,linestyle=":")
    scaled.text(common_end*0.96,86,"common data\nrange ends",ha="right",va="top",fontsize=8,color="#555555")
    rate.set(title="C  Rate per gram at equal progress",xlabel="Target acid conversion (%)",
             ylabel="Average specific rate / 10% loading",xlim=(8,42),ylim=(0.68,1.82))
    rate.axhline(1,color="#777777",linewidth=1.1,linestyle="--")
    rate.set_xticks([10,20,30,40])
    for ax in (actual,scaled,rate):
        ax.grid(alpha=.20)
        ax.spines[["top","right"]].set_visible(False)
    figure.suptitle("Pentyl propanoate at 40 °C, nominal alcohol:acid = 2:1",fontsize=14,weight="bold")
    figure.savefig(path,facecolor="white")


def plot_extrapolation(summary, path):
    case = read_json(summary)
    transfer = case["bed_transfer"]
    ref = transfer["batch_reference_concentration_g_L"]
    lo = min(row["nominal_catalyst_concentration_g_L"] for row in case["batch_condition_summary"].values())
    hi = transfer["largest_batch_concentration_g_L"]
    bulk = transfer["bed_concentrations_g_L"]["total_volume_proxy"]
    liquid = transfer["bed_concentrations_g_L"]["void_liquid_basis"]
    x = np.geomspace(.85*lo,1.10*liquid,301)
    u = x/ref
    all_fit = case["batch_fit"]["variants"]["printed"]["full_fit"]
    high_fit = case["batch_fit"]["high_pair_fit_lower_loadings_held_out"]["cases"]
    diagnostic = case["batch_fit"]["diagnostic_individual_loading_at_common_published_equilibrium"]
    # The plotted curves are explicit evaluations of the reported scalar
    # rate laws. They do not solve or refit the reaction model.
    fig = Figure(figsize=(9.1,5.2),layout="constrained",dpi=160)
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    ax.axvspan(lo,hi,color="#e2edf6",alpha=.9,label="Measured loading range")
    ax.axhline(1,color="#444444",linestyle=":",linewidth=1.2,label="Proportional / full-fit crowding")
    ax.plot(x,u**(all_fit["power"]["shape"]-1),color="#bb5e19",linewidth=2,label="Power: all four loadings")
    ax.plot(x,u**(high_fit["power"]["shape"]-1),color="#2478b4",linestyle="--",linewidth=2,
            label="Power: highest two only")
    ax.plot(x,1/(1+high_fit["crowding"]["shape"]*(u-1)),color="#279265",linestyle="--",linewidth=2,
            label="Crowding: highest two only")
    base_rate = diagnostic["10"]["k_ref_mol_kg_s"]
    conditions = case["batch_condition_summary"]
    ax.scatter([row["nominal_catalyst_concentration_g_L"] for row in conditions.values()],
               [diagnostic[name]["k_ref_mol_kg_s"]/base_rate for name in conditions],
               color="#171717",marker="o",zorder=5,s=36,label="Individual fits, common published equilibrium")
    for xpos, label in ((bulk,"Bed / total volume"),(liquid,"Bed / void liquid")):
        ax.axvline(xpos,color="#666666",linewidth=1.2,linestyle="-.")
        ax.text(xpos/1.015,1.48,label,rotation=90,va="top",ha="right",fontsize=9,color="#555555",
                bbox=dict(facecolor="white",alpha=.85,edgecolor="none",pad=1))
    ax.set(xscale="log",xlim=(x[0],x[-1]),ylim=(.55,1.51),
           xlabel="Catalyst concentration (g/L)",ylabel="Rate per gram / rate at 10% loading",
           title="Loading-law extrapolation to the packed bed is not identified")
    ax.set_xticks([lo,ref,hi,100,bulk,liquid],labels=[f"{lo:.1f}",f"{ref:.1f}",f"{hi:.1f}","100",f"{bulk:.0f}",f"{liquid:.0f}"])
    ax.grid(axis="y",alpha=.2)
    ax.legend(loc="upper left",fontsize=8,framealpha=.92)
    ax.spines[["top","right"]].set_visible(False)
    fig.savefig(path,facecolor="white")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--summary",type=Path)
    parser.add_argument("--extrapolation-output",type=Path)
    args = parser.parse_args()
    if bool(args.summary) != bool(args.extrapolation_output):
        parser.error("Supply both --summary and --extrapolation-output for the loading-law figure.")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    plot(args.output)
    if args.summary:
        args.extrapolation_output.parent.mkdir(parents=True,exist_ok=True)
        plot_extrapolation(args.summary,args.extrapolation_output)
