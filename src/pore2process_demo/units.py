"""Display conversions; the bridge and saved numerical values always use SI."""
from .bridge import NetworkInputs


def from_display(nx=6, ny=4, nz=4, spacing_mm=1., radius_um=50., viscosity_mPa_s=1., p_in_Pa=100., p_out_Pa=0.):
    return NetworkInputs(nx, ny, nz, spacing_mm*1e-3, radius_um*1e-6, viscosity_mPa_s*1e-3, p_in_Pa, p_out_Pa)


def microlitres_per_minute(flow_m3_s):
    return None if flow_m3_s is None else flow_m3_s * 60e9


def millilitres_per_minute(flow_m3_s):
    return None if flow_m3_s is None else flow_m3_s * 60e6
