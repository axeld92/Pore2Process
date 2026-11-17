"""Transfer the existing carbon-access solver to published bed lengths."""
from functools import lru_cache
from pathlib import Path

import numpy as np

from examples.adsorption.carbon_access import kernel as access
from examples.adsorption.disordered_memory_column import run as column
from pore2process_demo.bridge import runtime


def source():
    text = access.source()

    def replace(before,after):
        nonlocal text
        if text.count(before)!=1:
            raise ValueError(f"Missing unique transfer boundary: {before}")
        text = text.replace(before,after)

    replace("module CarbonAccessColumn","module CarbonAccessTransfer")
    replace("seam_history=false,direct_fraction", "seam_history=false,terminal_sections=0,direct_fraction")
    replace("    tau = sum(dot(r.share,r.seconds) for r in rules)\n    scale = liquid_volume_m3/(Q*blocks*tau)",
        """    0 <= terminal_sections < nsample || throw(ArgumentError("Invalid terminal section count"))
    kinds = [repeat(1:nsample,blocks); 1:terminal_sections]
    tau = sum(dot(rules[k].share,rules[k].seconds) for k in kinds)
    scale = liquid_volume_m3/(Q*tau)""")
    replace("    ns = nsample*blocks\n    kinds = repeat(1:nsample,blocks)", "    ns = length(kinds)")
    replace("out[step+1] = incoming[1]", "out[step+1] = dot(flow,c)")
    replace("liquid_volume_m3=Q*blocks*tau*scale", "liquid_volume_m3=Q*tau*scale")
    return text


@lru_cache(maxsize=1)
def solver():
    access.load()
    jl = runtime()
    jl.include_string(jl,source(),"CarbonAccessTransfer.jl")
    return jl.CarbonAccessTransfer.solve_column


@lru_cache(maxsize=4)
def transport(flow,probe_path):
    if flow==43.8:
        return column.transport("d800_primary",4)
    base = Path(probe_path)/f"q{round(10*flow)}"
    loader,_ = column.julia_models()
    paths = [str(base/name) for name in ("training.f64","operator.f64","rules.f64")]
    return loader(*paths,4)


def solve(flow,probe_path,p,scenario):
    model = transport(flow,probe_path)
    r = solver()(model.sections,model.flow,model.area,**p,**scenario)
    jl = runtime()
    return {str(f):np.array(v,dtype=float).tolist() if bool(jl.isa(v,jl.AbstractArray)) else float(v)
            for f in jl.fieldnames(jl.typeof(r)) for v in (getattr(r,str(f)),)}
