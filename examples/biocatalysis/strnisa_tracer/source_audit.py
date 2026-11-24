"""Check what absolute-time scale is represented by the normalized Fig. 3a.

This is a reading audit, not a claim about unpublished experiment processing.
"""
import tomllib
import numpy as np
from . import loop, native
from pore2process_demo.io import sha

FIG3A=loop.MARKERS.with_name("rtd_fig3a_audit.toml")


def audit():
    d=tomllib.loads(FIG3A.read_text())
    a=np.array(d["markers"],dtype=float)
    times=a[:,0]*5.
    x=(a[:,1]-d["axis_time_0_px"])/(d["axis_time_1_px"]-d["axis_time_0_px"])
    y=(d["axis_conc_0_px"]-a[:,2])/(d["axis_conc_0_px"]-d["axis_conc_1_px"])
    raw_t,raw_c,_=loop.digitized("A")
    yb=np.interp(times,raw_t,raw_c)/raw_c.max()
    inferred=float(np.dot(times,times)/np.dot(times,x))
    measured=loop.moments(*loop.measured("A"))["mean_s"]
    bed=60*next(row[6] for row in native.source()["table1"]["rows"] if row[0]=="A")
    return dict(fig3a_sha256=sha(FIG3A),fig3b_sha256=sha(loop.MARKERS),
                markers=len(a),inferred_time_scale_s=inferred,
                full_measured_loop_mean_s=measured,printed_A_reactor_mean_s=bed,
                normalized_time_max_error_using_full_loop_mean=float(np.max(np.abs(x-times/measured))),
                normalized_time_max_error_using_printed_bed_mean=float(np.max(np.abs(x-times/bed))),
                normalized_concentration_max_error=float(np.max(np.abs(y-yb))),
                normalized_concentration_rms_error=float(np.sqrt(np.mean((y-yb)**2))),
                inference="Fig. 3a experimental intensity sequence is consistent with Fig. 3b A divided by its peak, and time scaled by the complete loop mean. Overlay of normalized profiles alone cannot validate reactor-only dispersion.")


if __name__=="__main__":
    import json
    print(json.dumps(audit(),indent=2))
