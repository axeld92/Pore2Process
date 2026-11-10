"""Record actual hybrid reference errors and cold/warm timings, without a UI."""
import json
from pathlib import Path
from dataclasses import replace

import numpy as np
from pore2process_demo.bridge import NetworkInputs, predict_column_flow, run_case, runtime_versions

first = run_case()
warm = run_case()
largest = run_case(NetworkInputs(nx=20,ny=8,nz=8))
expected_Q = np.pi*(50e-6)**4*16*100/(8*1e-3*.005)
record = {
    "versions": runtime_versions(),
    "default": {
        "Q_in_m3_s": first.Q_in_m3_s, "Q_out_m3_s": first.Q_out_m3_s,
        "permeability_m2": first.permeability_m2, "column_flow_m3_s": predict_column_flow(first),
        "relative_flow_error": abs(first.Q_in_m3_s/expected_Q-1),
        "max_pressure_error_Pa": float(np.max(np.abs(first.pressures_Pa-(100-100*first.coordinates_m[:,0]/.005)))),
        "max_interior_imbalance_m3_s": first.max_interior_imbalance_m3_s,
        "flow_mismatch_m3_s": first.flow_mismatch_m3_s,
        "normalized_conservation_error": first.normalized_conservation_error,
    },
    "timings_s": {"cold_runtime_initialization": first.runtime_initialization_s,
                  "first_bridge_call_including_compilation": first.bridge_call_s,
                  "first_julia_core": first.total_s,
                  "warm_default_bridge_call": warm.bridge_call_s, "warm_default_julia_core": warm.total_s,
                  "largest_ui_bridge_call": largest.bridge_call_s, "largest_ui_julia_core": largest.total_s},
    "largest_ui_counts": {"nodes": len(largest.pressures_Pa), "throats": len(largest.edges)},
    "note": "Single local observations after dependency precompilation, not a performance benchmark; excludes graphics/export time.",
}
destination = Path('examples/hydraulics/capillary_network/reference/verified_reference.json')
destination.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
print(destination.read_text())
