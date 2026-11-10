"""Inspect actual native exports and call the shared Julia numerical package."""
import json
from pathlib import Path

import numpy as np
from pore2process_demo.io import write_json

from pore2process_demo.paths import ROOT, project_path
def read_protocol(path=None):
    return json.loads((project_path(Path(path).absolute()) if path else ROOT/"examples/filtration/graded_support/protocol.json").read_text())

def profile(protocol,bulk):
    g=protocol["geometry"]
    x=np.array(g["body_knots_fraction"])*g["body_length_m"]
    r=np.array([bulk if v=="bulk" else v*g["neck_radius_m"] for v in g["radius_factors"]])
    return x,r

def load_case(jl,R,R0,protocol,Z):
    p=protocol["process"]
    v=jl.graded_support_load(R,R0,p["baseline_support_flux_at_1bar_m_h"]/3600/1e5,
        Z,p["filtration_flux_LMH"]/3.6e6,p["membrane_area_m2"],p["pump_efficiency_fraction"])
    return {name:float(getattr(v,name)) for name in (
        "measured_body_resistance_ratio","support_permeance_m_Pa_s","support_impedance_Pa_s_m",
        "non_support_impedance_Pa_s_m","pressure_drop_Pa","flow_m3_s","shaft_power_W",
        "specific_shaft_energy_kWh_m3","non_support_fraction","support_fraction")}

def analyze(directory,raw,jl,protocol):
    directory=Path(directory)
    planes=np.atleast_1d(np.genfromtxt(directory/"planes.csv",delimiter=",",names=True))
    samples=np.atleast_1d(np.genfromtxt(directory/"midplane.csv",delimiter=",",names=True,dtype=None,encoding="utf-8"))
    history=np.atleast_1d(np.genfromtxt(directory/"history.csv",delimiter=",",names=True))
    for data in (planes,samples,history):
        if not all(np.all(np.isfinite(data[name])) for name in data.dtype.names):
            raise RuntimeError("Non-finite native CSV fields.")
    if len(planes)!=raw["axial_intervals"]+1 or not np.array_equal(planes["i"],np.arange(len(planes))):
        raise RuntimeError("Axial plane shape/order differs from native contract.")
    dx=raw["dx_m"]
    if not np.allclose(planes["x_m"],planes["i"]*dx,rtol=0,atol=1e-12*dx):
        raise RuntimeError("Plane coordinates disagree with zero-based indices.")
    for key,coord in (("i","x_m"),("j","y_m"),("k","z_m")):
        if samples[key].dtype!=np.dtype("int64"):
            raise RuntimeError("Native lattice indices must materialize as int64.")
        shift=0 if key=="i" else .5
        if not np.allclose(samples[coord],(samples[key]+shift)*dx,rtol=0,atol=1e-11*dx):
            raise RuntimeError("Native coordinate convention changed.")
    middle=raw["axial_intervals"]//2
    if not np.all(samples["i"]==middle) or len(samples)!=int(planes["node_count"][middle]):
        raise RuntimeError("Midplane shapes/indices differ from the section geometry.")
    i1=round(raw["lead_length_each_m"]/dx)
    i2=round((raw["lead_length_each_m"]+raw["body_length_m"])/dx)
    dp=float(planes["mean_pressure_Pa"][i1]-planes["mean_pressure_Pa"][i2])
    q=float(planes["Q_m3_s"][middle])
    for a,b in ((dp,raw["tap_delta_p_Pa"]),(q,raw["Q_mid_m3_s"])):
        if not np.isclose(a,b,rtol=1e-13,atol=1e-30):
            raise RuntimeError("Native JSON and plane CSV measurements differ.")
    if dp*q<=0 or abs(raw["inlet_pressure_Pa"]-raw["imposed_delta_p_Pa"])>1e-8 or abs(raw["outlet_pressure_Pa"])>1e-8:
        raise RuntimeError("Pressure/flow sign or boundary snapshot differs.")
    A=raw["cell_width_m"]**2; mu=raw["viscosity_Pa_s"]; L=raw["body_length_m"]
    k=float(jl.infer_permeability(L,A,mu,dp,q))
    R=float(jl.predict_column_pressure_drop(k,mu,L,A,1.))
    x,r=profile(protocol,raw["bulk_radius_m"])
    vector=jl.seval("Vector{Float64}")
    reference=jl.graded_pore_reference(vector(x.tolist()),vector(r.tolist()),mu,A)
    mass=planes["mass_flow_kg_s"]
    scale=float(np.mean(np.abs(mass)))
    speed=float(np.max(planes["max_speed_m_s"]))
    m={"permeability_cell_m2":k,"body_resistance_Pa_s_m3":R,"body_impedance_Pa_s_m":A*R,
       "body_pressure_drop_Pa":dp,"measured_flow_m3_s":q,
       "all_plane_mass_relative_spread":float(np.ptp(mass)/scale),
       "inlet_outlet_mass_relative_difference":float(abs(mass[0]-mass[-1])/scale),
       "maximum_speed_m_s":speed,"max_Re":speed*raw["diameter_m"]*raw["reference_density_kg_m3"]/mu,
       "max_Mach":float(speed*raw["dt_s"]/dx*np.sqrt(3)),
       "max_density_departure":float(np.max(planes["max_density_departure"])),
       "lubrication_reference_resistance_Pa_s_m3":float(reference.hydraulic_resistance_Pa_s_m3),
       "lubrication_relative_difference":float(R/reference.hydraulic_resistance_Pa_s_m3-1),
       "continuous_body_volume_porosity":float(reference.volume_porosity)}
    if raw["bulk_radius_m"]==protocol["geometry"]["neck_radius_m"]:
        r2=(samples["y_m"]-raw["mesh_span_m"]/2)**2+(samples["z_m"]-raw["mesh_span_m"]/2)**2
        exact=dp/L/(4*mu)*((raw["diameter_m"]/2)**2-r2)
        m["baseline_profile_relative_L2"]=float(np.linalg.norm(samples["ux_m_s"]-exact)/np.linalg.norm(exact))
    for field,key in (("Q_mid_m3_s","flow"),("delta_p_taps_Pa","pressure")):
        window=history[field][-raw["convergence_window_samples"]:]
        if len(window)!=raw["convergence_window_samples"]:
            raise RuntimeError("Native convergence window is incomplete.")
        value=float(np.std(window,ddof=1)/abs(np.mean(window)))
        m[f"convergence_{key}_relative_std"]=value
        if value>raw["steady_tolerance"]*(1+1e-5):
            raise RuntimeError("Exported history fails the stock OpenLB ValueTracer criterion.")
    if not all(np.isfinite(list(m.values()))):
        raise RuntimeError("Non-finite analysis values.")
    limits=protocol["acceptance"]
    checks={name:m[name]<=limits[name] for name in (
        "all_plane_mass_relative_spread","inlet_outlet_mass_relative_difference")}
    checks["low_Re"]=m["max_Re"]<limits["native_max_Re"]
    checks["low_Mach"]=m["max_Mach"]<limits["native_max_Mach"]
    if "baseline_profile_relative_L2" in m:
        checks["analytical_k"]=abs(m["lubrication_relative_difference"])<limits["baseline_analytical_k_relative_error"]
        checks["analytical_profile"]=m["baseline_profile_relative_L2"]<limits["baseline_velocity_profile_relative_L2"]
    if not all(checks.values()):
        raise RuntimeError(f"Native acceptance failed at {directory}: {checks}; metrics={m}")
    return m,{name:bool(v) for name,v in checks.items()}
