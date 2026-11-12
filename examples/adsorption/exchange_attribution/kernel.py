"""Extend the pinned exchange solver by separating carbon and liquid shares."""
from functools import lru_cache

import numpy as np

from examples.adsorption.exchange_column.run import JULIA, solver as baseline_solver
from pore2process_demo.bridge import runtime
from pore2process_demo.io import sha

ORIGINAL = JULIA / "ExchangeColumn.jl"
ORIGINAL_SHA = "614863a566edf41d432a26a54cf3f3ba3f86a1e1c04823a357e00267eee23a52"


def source():
    if sha(ORIGINAL)!=ORIGINAL_SHA:
        raise ValueError("Pinned exchange-column solver changed")
    code = ORIGINAL.read_text()
    changes = (
        ("module ExchangeColumn\n", "module AllocationColumn\n"),
        ("startup_dt_s=.05,kwargs...)", "startup_dt_s=.05,carbon_fractions=nothing,kwargs...)"),
        ("    M = solid.*share",
         """    if carbon_fractions === nothing
        carbon_share = share
    else
        carbon = Vector{Float64}(carbon_fractions)
        length(carbon)==n && all(isfinite,carbon) && all(>=(0),carbon) &&
            abs(sum(carbon)-1)<=1e-12 || throw(ArgumentError("Carbon shares must be nonnegative and sum to one"))
        carbon_share = carbon/(slices*section_count)
    end
    M = solid.*carbon_share"""),
        ("      liquid_volume_m3=liquid,physical_adsorbent_mass_kg=solid,",
         """      liquid_volume_m3=liquid,physical_adsorbent_mass_kg=solid,
      class_carbon_mass_kg=M*(slices*section_count),class_liquid_volume_m3=V*(slices*section_count),"""),
    )
    for old, new in changes:
        if code.count(old)!=1:
            raise ValueError("Pinned solver has an unexpected edit boundary")
        code = code.replace(old,new)
    return code


@lru_cache(maxsize=1)
def load(path):
    baseline_solver()
    jl = runtime()
    jl.include(str(path))
    return jl.AllocationColumn.solve_column


def solve(path, shares, delays, mixing, carbon, **p):
    jl = runtime()
    r = load(path)(np.array(shares),np.array(delays),np.array(mixing),carbon_fractions=np.array(carbon),**p)
    out = {}
    for field in jl.fieldnames(jl.typeof(r)):
        key = str(field)
        value = getattr(r,key)
        out[key] = np.array(value,dtype=np.float64,copy=True).tolist() if bool(jl.isa(value,jl.AbstractArray)) else float(value)
    return out


def verify(path):
    p = dict(periods=2,subcells=2,length_m=.02,area_m2=1e-4,porosity=.4,
             superficial_velocity_m_s=.02,inlet_concentration_kg_m3=.15,
             particle_diameter_m=.0008,particle_density_kg_m3=1290.,film_m_s=4.5e-5,
             diffusivity_m2_s=1.85e-12,radial_cells=12,t_end_s=30.,dt_s=.1)
    w = np.array([.5,.5])
    delays = np.array([.6,1.4])
    mix = np.array([[.8,.2],[.2,.8]])
    ref = baseline_solver()(w,delays,mix,**p)
    retained = solve(path,w,delays,mix,w*delays,**p)
    keys = ("outlet_concentration_kg_m3","adsorbed_inventory_kg","liquid_inventory_kg")
    for name in keys:
        np.testing.assert_array_equal(retained[name],np.array(getattr(ref,name)))
    moved = solve(path,w,delays,mix,w,**p)
    np.testing.assert_allclose(np.array(moved["class_carbon_mass_kg"])/moved["physical_adsorbent_mass_kg"],w)
    np.testing.assert_array_equal(retained["class_liquid_volume_m3"],moved["class_liquid_volume_m3"])
    if np.array_equal(retained[keys[0]],moved[keys[0]]):
        raise ValueError("Redistribution did not affect adsorption")
    p["diffusivity_m2_s"] = 0.
    tracer_a = solve(path,w,delays,mix,w*delays,**p)
    tracer_b = solve(path,w,delays,mix,w,**p)
    np.testing.assert_array_equal(tracer_a[keys[0]],tracer_b[keys[0]])
    return dict(baseline_reproduced=True, carbon_moves_without_liquid=True,
                adsorption_response_changes=True, inert_response_preserved=True)
