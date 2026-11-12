"""Preserve the prior shared-carbon solver while adding per-path film coefficients."""
from functools import lru_cache

import numpy as np

from examples.adsorption.shared_contact.run import JULIA, solver as shared_solver
from pore2process_demo.bridge import runtime
from pore2process_demo.io import sha

ORIGINAL = JULIA/"SharedContactColumn.jl"
ORIGINAL_SHA = "d64b4554de5c42e544534b17317eb2a78777edbf9884f2d2866860afe217b2e3"


def source():
    if sha(ORIGINAL)!=ORIGINAL_SHA:
        raise ValueError("The archived shared-carbon solver changed")
    code = ORIGINAL.read_text()
    edits = (
        ("module SharedContactColumn\n","module FilmByPathColumn\n"),
        ("linear_partition_m3_kg=nothing,kwargs...)","linear_partition_m3_kg=nothing,film_by_path_m_s=nothing,kwargs...)"),
        ("            cs[i] = surface_value", """            if beta[i]==0
                # A zero film coefficient makes this surface sector reflective.
                qs[i] = q
                mean += weights[i]*q
                derivative += weights[i]
                continue
            end
            cs[i] = surface_value"""),
        ("    active = p.kf>0 && p.Ds>0 && cap>0\n    conductance = active ? area*p.kf.*weights : zeros(n)",
         """    film_by_path = film_by_path_m_s===nothing ? fill(p.kf,n) : Vector{Float64}(film_by_path_m_s)
    length(film_by_path)==n && all(isfinite,film_by_path) && all(>=(0),film_by_path) ||
        throw(ArgumentError("Path film coefficients must be finite and nonnegative"))
    film_scale = maximum(film_by_path)
    active = any(>(0),film_by_path.*weights) && p.Ds>0 && cap>0
    conductance = active ? area.*film_by_path.*weights : zeros(n)"""),
        ("FilmParticle.operator(p.nr,p.R,p.Ds,p.kf,h)",
         "FilmParticle.operator(p.nr,p.R,p.Ds,active ? film_scale : 0.,h)"),
        ("        beta = active ? (ks/p.kf).*denom./liquid_denom : zeros(n)",
         """        beta = [active && film_by_path[i]>0 && weights[i]>0 ?
                ks/film_by_path[i]*denom[i]/liquid_denom[i] : 0. for i in 1:n]"""),
        ("if active && weights[i]>0\n                        concentration", "if active && conductance[i]>0\n                        concentration"),
        ("film = p.kf*(concentration-cs[i])", "film = film_by_path[i]*(concentration-cs[i])"),
        ("normalized_interface_flux_mismatch=maxflux/max(p.kf*cin,1e-30)",
         "normalized_interface_flux_mismatch=maxflux/max(film_scale*cin,1e-30)"),
        ("      total_film_conductance_m3_s=ns*sum(conductance),",
         """      total_film_conductance_m3_s=ns*sum(conductance),
      film_coefficient_by_path_m_s=film_by_path,film_conductance_by_path_m3_s=ns.*conductance,"""),
    )
    for old,new in edits:
        if code.count(old)!=1:
            raise ValueError(f"Expected one archived edit boundary: {old}")
        code=code.replace(old,new)
    return code


@lru_cache(maxsize=1)
def load(path):
    shared_solver()
    jl = runtime()
    jl.include(str(path))
    return jl.FilmByPathColumn.solve_column


def solve(path,cl,mix,weights,film,p):
    jl = runtime()
    result = load(path)(np.array(cl["flow_fractions"]),np.array(cl["normalized_transit_times"]),
                        np.array(mix),np.array(weights),film_by_path_m_s=np.array(film),**p)
    out = {}
    for field in jl.fieldnames(jl.typeof(result)):
        key = str(field)
        value = getattr(result,key)
        out[key] = np.array(value,dtype=float,copy=True).tolist() if bool(jl.isa(value,jl.AbstractArray)) else float(value)
    return out


def verify(path):
    cl = dict(flow_fractions=[.5,.5],normalized_transit_times=[.5,1.5])
    weights = [.3,.7]
    mix = np.array([[.8,.2],[.2,.8]])
    p = dict(periods=2,subcells=2,length_m=.02,area_m2=1e-4,porosity=.4,
             superficial_velocity_m_s=.02,inlet_concentration_kg_m3=.15,
             particle_diameter_m=.0008,particle_density_kg_m3=1290.,film_m_s=4.5e-5,
             diffusivity_m2_s=1.85e-12,radial_cells=12,t_end_s=30.,dt_s=.1)
    old = shared_solver()(np.array(cl["flow_fractions"]),np.array(cl["normalized_transit_times"]),
                         mix,np.array(weights),**p)
    baseline = solve(path,cl,mix,weights,[p["film_m_s"]]*2,p)
    for name in ("outlet_concentration_kg_m3","adsorbed_inventory_kg","liquid_inventory_kg"):
        np.testing.assert_allclose(baseline[name],np.array(getattr(old,name)),rtol=1e-10,atol=1e-15)
    varied = solve(path,cl,mix,weights,[0.,8e-5],p)
    assert varied["normalized_mass_residual"]<1e-9 and varied["normalized_interface_flux_mismatch"]<1e-9
    assert varied["film_conductance_by_path_m3_s"][0]==0. and varied["contact_area_by_path_m2"][0]>0.
    assert varied["adsorbed_inventory_kg"][-1]>0.
    assert not np.array_equal(varied["adsorbed_inventory_kg"],baseline["adsorbed_inventory_kg"])
    stopped = solve(path,cl,mix,weights,[0.,0.],p)
    inert = solve(path,cl,mix,weights,[2e-5,8e-5],{**p,"diffusivity_m2_s":0.})
    np.testing.assert_array_equal(stopped["outlet_concentration_kg_m3"],inert["outlet_concentration_kg_m3"])
    return dict(same_film_recovers_parent=True,zero_film_with_positive_area_reflects=True,
                nonuniform_transfer_conservative=True,inert_carrier_unchanged=True)
