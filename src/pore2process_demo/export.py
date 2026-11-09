"""Reproducible SI-valued exports from one immutable completed-run snapshot."""
from dataclasses import asdict
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from .bridge import ColumnInputs, predict_column_flow, runtime_versions
from .visualization import network_mesh, pressure_figure


def summary(result, column=ColumnInputs()):
    return {
        "schema_version": 1, "run_id": result.run_id, "timestamp_utc": result.timestamp_utc,
        "model": "Steady idealized cylindrical-capillary lattice; uniform viscosity; no gravity or junction loss",
        "inputs": asdict(result.inputs),
        "units": {"coordinates": "m", "spacing_m": "m", "radius_m": "m", "viscosity_Pa_s": "Pa s",
                  "pressure": "Pa", "flow": "m^3/s", "permeability": "m^2", "area": "m^2", "timing": "s"},
        "indexing": "Zero-based node_id = i + nx*(j + ny*k); x varies fastest",
        "flow_orientation": "Positive from first endpoint to second endpoint; edges emitted positive x, y, z per node",
        "domain": {"length_m": result.length_m, "representative_area_m2": result.area_m2,
                   "y_bounds_m": [0,result.inputs.ny*result.inputs.spacing_m],
                   "z_bounds_m": [0,result.inputs.nz*result.inputs.spacing_m],
                   "area_convention": "ny*nz*spacing_m^2; one transverse cell per node"},
        "counts": {"nodes": len(result.pressures_Pa), "throats": len(result.flows_m3_s)},
        "Q_in_m3_s": result.Q_in_m3_s, "Q_out_m3_s": result.Q_out_m3_s,
        "permeability_m2": result.permeability_m2,
        "column": {"network_run_id": result.run_id, "inputs": asdict(column),
                   "viscosity_Pa_s": result.inputs.viscosity_Pa_s, "flow_m3_s": predict_column_flow(result,column),
                   "assumption": "Same idealized repeating medium at column scale; demonstration of upscaling"},
        "diagnostics": {name: getattr(result,name) for name in ("max_interior_imbalance_m3_s", "flow_mismatch_m3_s", "normalized_conservation_error")},
        "timing_s": {"julia_runtime_initialization": result.runtime_initialization_s,
                     "bridge_call_including_first_compilation_and_conversion": result.bridge_call_s,
                     "julia_core": result.total_s, "linear_solve": result.linear_solve_s},
        "versions": runtime_versions(),
    }


def write_exports(result, output, column=ColumnInputs()):
    # Validate and build summary before writing any files.
    rendered = json.dumps(summary(result,column), indent=2, allow_nan=False)
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True)
    (output / "summary.json").write_text(rendered + "\n")
    nodes = np.column_stack((np.arange(len(result.pressures_Pa)), result.coordinates_m, result.pressures_Pa))
    np.savetxt(output / "nodes.csv", nodes, delimiter=",", fmt=["%d"]+["%.17g"]*4,
               header="node_id_zero_based,x_m,y_m,z_m,pressure_Pa", comments="")
    throats = np.column_stack((result.edges,result.radii_m,result.lengths_m,result.flows_m3_s))
    np.savetxt(output / "throats.csv", throats, delimiter=",", fmt=["%d","%d"]+["%.17g"]*3,
               header="from_node_id_zero_based,to_node_id_zero_based,radius_m,length_m,flow_from_to_m3_s", comments="")
    figure = pressure_figure(result)
    figure.savefig(output / "pressure_profile.png", dpi=160)
    figure.clear()
    network_mesh(result).save(output / "network.vtp")
    return output


def zip_export(result, column=ColumnInputs()):
    """Private temporary files, returned as an in-memory session-specific ZIP."""
    buffer = io.BytesIO()
    with TemporaryDirectory(prefix="pore2process-") as temporary:
        directory = write_exports(result,temporary,column)
        with ZipFile(buffer,"w",ZIP_DEFLATED) as archive:
            for path in sorted(directory.iterdir()):
                archive.write(path,path.name)
    buffer.seek(0)
    return buffer
