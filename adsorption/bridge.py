"""Load the new Julia module without changing the frozen filtration package."""
from functools import lru_cache
import numpy as np
from pore2process_demo.bridge import ROOT,runtime,_check_thread

@lru_cache(maxsize=1)
def _loaded_module():
    jl=runtime()
    jl.include(str(ROOT/'julia/adsorption/AdsorptionColumn.jl'))
    return jl.AdsorptionColumn

def module():
    _check_thread()
    return _loaded_module()

def column(**inputs):
    _check_thread()
    jl=runtime();v=module().solve_column(**inputs);out={}
    for field in jl.fieldnames(jl.typeof(v)):
        k=str(field);x=getattr(v,k)
        if x is None:out[k]=None
        elif k in ('x_m','time_s','outlet_concentration_kg_m3','final_concentration_kg_m3','final_loading_kg_kg',
                'inventory_kg','adsorbed_inventory_kg','liquid_inventory_kg','cumulative_inlet_kg','cumulative_outlet_kg','mass_residual_kg'):
            out[k]=np.array(x,dtype=np.float64,copy=True).tolist()
        elif k in ('cells','steps','max_Newton_iterations'):out[k]=int(x)
        else:out[k]=float(x)
    return out
