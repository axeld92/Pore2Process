"""Scientific plots and SI wall meshes from the saved native study."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from examples.filtration.graded_support.model import profile

def render(path):
    import pyvista as pv
    path=Path(path); base=path.parent
    s=json.loads(path.read_text()); p=s["protocol"]
    colors=("#275d8a","#c36523","#32816b")
    names=("Straight, 4 µm","Interior 5 µm","Interior 6 µm")
    fig=plt.figure(figsize=(11.8,8.5),layout="constrained")
    theta=np.linspace(0,2*np.pi,65)
    for row,(support,label,color) in enumerate(zip(s["supports"],names,colors,strict=True)):
        x,r=profile(p,support["bulk_radius_m"])
        ax=fig.add_subplot(3,2,2*row+1)
        ax.fill_between(x*1e6,-r*1e6,r*1e6,color=color,alpha=.3)
        ax.plot(x*1e6,r*1e6,c=color); ax.plot(x*1e6,-r*1e6,c=color)
        ax.axvspan(32,40,color="#666666",alpha=.12)
        ax.text(36,6.7,"Fixed outlet\nneck / interface",ha="center",va="top",fontsize=7)
        ax.set(xlim=(0,40),ylim=(-7,7),xlabel="Support x (µm)",ylabel="Radial position (µm)",title=label)
        ax.set_aspect("equal")
        raw=support["grids"][-1]["raw"]
        planes=np.genfromtxt(base/support["grids"][-1]["case_id"]/"planes.csv",delimiter=",",names=True)
        xx=np.linspace(0,x[-1],121)
        rr=np.interp(xx,x,r)
        X=np.repeat(xx[:,None],len(theta),axis=1)
        Y=p["geometry"]["representative_cell_width_m"]/2+rr[:,None]*np.cos(theta)
        Z=p["geometry"]["representative_cell_width_m"]/2+rr[:,None]*np.sin(theta)
        pressure=np.interp(xx+raw["lead_length_each_m"],planes["x_m"],planes["mean_pressure_Pa"])
        # This is a visual wall surface; pressure is the section mean, not a
        # computed wall-pressure or wall-stress field. Native coordinates use
        # span/2; shift to the centre of the declared representative cell.
        grid=pv.StructuredGrid(X,Y,Z)
        grid["section_mean_pressure_Pa"]=np.broadcast_to(pressure[:,None],X.shape).ravel(order="F")
        grid.field_data["metadata"]=["SI wall visualization; section-mean pressure interpolated from OpenLB planes; numerical leads omitted; no wall stress inferred"]
        grid.field_data["native_to_visual_offset_m"]=[-raw["lead_length_each_m"],4e-6,4e-6]
        grid.extract_surface(algorithm="dataset_surface").save(base/f"{support['name']}_wall.vtp")
        ax=fig.add_subplot(3,2,2*row+2,projection="3d")
        faces=plt.cm.viridis(np.broadcast_to((pressure/10)[:,None],X.shape))
        ax.plot_surface(X*1e6,(Y-1e-5)*1e6,(Z-1e-5)*1e6,facecolors=faces,shade=False,linewidth=0,rstride=3,cstride=3)
        ax.set(xlabel="x (µm)",ylabel="y (µm)",zlabel="z (µm)",xlim=(0,40),ylim=(-7,7),zlim=(-7,7),
               xticks=[0,20,40],yticks=[-6,0,6],zticks=[-6,0,6])
        ax.set_box_aspect((40,14,14)); ax.view_init(elev=18,azim=-65)
        ax.set_title("Wall shape; colour is section-mean pressure",fontsize=9)
    scalar=plt.cm.ScalarMappable(norm=plt.Normalize(0,10),cmap="viridis")
    fig.colorbar(scalar,ax=fig.axes[1::2],label="Section-mean gauge pressure (Pa)",shrink=.75,pad=.08)
    fig.savefig(base/"geometries.png",dpi=170); plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11.5,8),layout="constrained")
    nominal=s["scenarios"][0]["cases"]
    low=s["scenarios"][1]["cases"]; high=s["scenarios"][2]["cases"]
    R0=s["supports"][0]["body_resistance_Pa_s_m3"]
    for support,label,color in zip(s["supports"],names,colors,strict=True):
        axes[0,0].plot([v["raw"]["resolution"] for v in support["grids"]],
            [v["metrics"]["body_resistance_Pa_s_m3"]/R0 for v in support["grids"]],"o-",c=color,label=label)
        fine=support["grids"][-1]; raw=fine["raw"]
        planes=np.genfromtxt(base/fine["case_id"]/"planes.csv",delimiter=",",names=True)
        mask=(planes["x_m"]>=raw["lead_length_each_m"]-1e-12)&(planes["x_m"]<=raw["lead_length_each_m"]+raw["body_length_m"]+1e-12)
        axes[0,1].plot((planes["x_m"][mask]-raw["lead_length_each_m"])*1e6,planes["mean_pressure_Pa"][mask],c=color,label=label)
    axes[0,0].set(xlabel="Nodes across the fixed neck diameter",ylabel="Native body resistance / fine straight baseline",title="Measured OpenLB resistance: three grids")
    axes[0,1].set(xlabel="Support x (µm)",ylabel="Section-mean gauge pressure (Pa)",title="Native pressure profiles; same imposed 10 Pa")
    values=np.array([c["actual_specific_shaft_energy_kWh_m3"] for c in nominal])
    bounds=np.array([[c["actual_specific_shaft_energy_kWh_m3"] for c in scenario] for scenario in (low,high)])
    axes[1,0].bar(names,values,color=colors)
    axes[1,0].errorbar(names,values,yerr=np.vstack((values-bounds.min(axis=0),bounds.max(axis=0)-values)),fmt="none",c="black",capsize=5,label="PDF reading sensitivity only")
    axes[1,0].set(ylabel="Actual DWSIM shaft energy (kWh/m³)",title="Literature-anchored scenario; fixed 100 LMH")
    shares=np.array([c["julia"]["non_support_fraction"] for c in nominal])
    axes[1,1].bar(names,100*(1-shares),label="Support",color="#4c86ab")
    axes[1,1].bar(names,100*shares,bottom=100*(1-shares),label="Effective non-support",color="#b7c8d5")
    axes[1,1].set(ylabel="Impedance share (%)",ylim=(0,105),title="Fixed effective non-support impedance")
    for ax in axes.flat:
        ax.grid(alpha=.2,axis="y"); ax.legend(fontsize=8)
    fig.savefig(base/"comparison.png",dpi=170); plt.close(fig)

if __name__=="__main__":
    import argparse
    a=argparse.ArgumentParser(description=__doc__); a.add_argument("summary",type=Path)
    render(a.parse_args().summary)
