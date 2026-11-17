"""Use the archived transport solver with a two-cohort particle boundary."""
from functools import lru_cache

import numpy as np

from examples.adsorption.disordered_memory_column import run as column
from pore2process_demo.bridge import runtime
from pore2process_demo.io import sha

ORIGINAL = column.HERE/"model.jl"
ORIGINAL_SHA = "4c16107a4d05baab77150fe82db4b676e339b2ab652c54d56269918409b1d4d2"


def source():
    if sha(ORIGINAL)!=ORIGINAL_SHA:
        raise ValueError("The archived transport solver changed")
    text = ORIGINAL.read_text()

    def replace(old,new):
        nonlocal text
        if text.count(old)!=1:
            raise ValueError(f"Expected one boundary in the archived solver: {old}")
        text = text.replace(old,new)

    def change_block(start,end,new):
        a = text.index(start)
        b = text.index(end,a)
        replace(text[a:b],new)

    replace("module MemoryColumn", "module CarbonAccessColumn")
    replace("using ..SharedContactColumn: particle_surface!", "using ..AccessCohorts")
    replace("seam_history=false,kwargs...)", "seam_history=false,direct_fraction=1.,slow_ratio=0.,kwargs...)")
    change_block("    cap = p.qe(cin)\n","    history = zeros", """    bed = AccessCohorts.setup(p,Matrix{Float64}(contact),ns,mass,area,Q.*flow,cin,
                              direct_fraction,slow_ratio)
    active = any(bed.active)
    c,incoming = zeros(n),zeros(n)
    solid_groups = zeros(2,nstep+1)
""")
    change_block("    ks = p.rho*p.Ds/op.half_width", "    t0 = time_ns()", "")
    replace("            copyto!(interior,particle)\n            ldiv!(op.F,interior)",
            "            AccessCohorts.diffuse!(bed,op)")
    change_block("            if active\n                @inbounds", "            @inbounds for i in 1:n\n                minc",
        """            if active
                r = AccessCohorts.contact!(bed,c,s,k,h,op,capacity_kg_kg,affinity_m3_kg,linear_partition_m3_kg)
                solid[step+1] += sum(r.adsorbed)
                solid_groups[1,step+1] += r.adsorbed[1]
                solid_groups[2,step+1] += r.adsorbed[2]
                maxq = max(maxq,r.maxq)
                maxflux = max(maxflux,r.maxflux)
                uptake += r.uptake
                release += r.release
            end
""")
    replace("adsorbed_inventory_kg=solid,mass_residual_kg=residual,",
        """adsorbed_inventory_kg=solid,mass_residual_kg=residual,
      adsorbed_inventory_by_cohort_kg=solid_groups,
      carbon_by_cohort_kg=dry_mass_kg.*collect(bed.fractions),
      equivalent_film_by_cohort_m_s=collect(bed.film),
      effective_conductance_m3_s=ns*area*sum(bed.film.*bed.fractions),
      equilibrium_adsorbed_ceiling_kg=dry_mass_kg*bed.cap*sum(bed.fractions.*bed.active),""")
    return text


@lru_cache(maxsize=1)
def load():
    column.julia_models()
    jl = runtime()
    jl.include(str(ORIGINAL.parent.parent/"carbon_access/cohorts.jl"))
    jl.include_string(jl,source(),"CarbonAccessColumn.jl")
    return jl.CarbonAccessColumn


def solve(size,p,scenario):
    model = column.transport(size+"_primary",4)
    r = load().solve_column(model.sections,model.flow,model.area,**p,**scenario)
    jl = runtime()
    return {str(f):np.array(v,dtype=float).tolist() if bool(jl.isa(v,jl.AbstractArray)) else float(v)
            for f in jl.fieldnames(jl.typeof(r)) for v in (getattr(r,str(f)),)}
