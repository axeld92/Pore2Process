"""The only JuliaCall configuration/import site. All public values use SI units."""
from dataclasses import asdict, dataclass, fields
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import version
import json
import numbers
import os
from pathlib import Path
import platform
import sys
import threading
from time import perf_counter
from uuid import uuid4

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def _check_thread():
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError("This prototype calls Julia only on the main server thread. Use its event-loop callbacks; background Julia calls are unsupported.")


def _input_error(error):
    # Keep actionable validation text in the product; logs retain traceback details.
    message = str(error).split("Stacktrace:",1)[0].strip()
    return message.removeprefix("ArgumentError: ")


@dataclass(frozen=True)
class NetworkInputs:
    nx: int = 6
    ny: int = 4
    nz: int = 4
    spacing_m: float = 1e-3
    radius_m: float = 50e-6
    viscosity_Pa_s: float = 1e-3
    p_in_Pa: float = 100.0
    p_out_Pa: float = 0.0


@dataclass(frozen=True)
class ColumnInputs:
    length_m: float = .02
    area_m2: float = 1e-4
    delta_p_Pa: float = 1000.0


@dataclass(frozen=True)
class CaseResult:
    inputs: NetworkInputs
    run_id: str
    timestamp_utc: str
    coordinates_m: np.ndarray
    edges: np.ndarray
    pressures_Pa: np.ndarray
    flows_m3_s: np.ndarray
    radii_m: np.ndarray
    lengths_m: np.ndarray
    conductances_m3_Pa_s: np.ndarray
    net_outflows_m3_s: np.ndarray
    inlet_nodes: np.ndarray
    outlet_nodes: np.ndarray
    interior_nodes: np.ndarray
    Q_in_m3_s: float
    Q_out_m3_s: float
    length_m: float
    area_m2: float
    permeability_m2: float | None
    max_interior_imbalance_m3_s: float
    flow_mismatch_m3_s: float
    normalized_conservation_error: float
    linear_solve_s: float
    total_s: float
    runtime_initialization_s: float
    bridge_call_s: float


@lru_cache(maxsize=1)
def runtime():
    """Start one persistent Julia runtime; never install dependencies on launch."""
    _check_thread()
    config_file = ROOT / ".local/julia_runtime.json"
    if not config_file.exists():
        raise RuntimeError("Julia bridge is not bootstrapped. Run: uv run python scripts/bootstrap_julia.py")
    config = json.loads(config_file.read_text())
    if config["pythoncall_version"] != version("juliacall"):
        raise RuntimeError("JuliaCall/PythonCall versions changed. Rerun scripts/bootstrap_julia.py.")
    if config["project"] != str(ROOT / "julia") or not Path(config["executable"]).is_file():
        raise RuntimeError("Project or Julia executable moved. Rerun scripts/bootstrap_julia.py.")
    if "juliacall" in sys.modules:
        raise RuntimeError("JuliaCall was imported before bridge configuration. Restart and import this bridge first.")
    for key, value in (("PYTHON_JULIACALL_EXE", config["executable"]), ("PYTHON_JULIACALL_PROJECT", str(ROOT / "julia"))):
        if key in os.environ and os.environ[key] != value:
            raise RuntimeError(f"{key} conflicts with the bootstrapped runtime. Unset it or rerun bootstrap.")
        os.environ[key] = value
    if config.get("depot_path"):
        os.environ["JULIA_DEPOT_PATH"] = config["depot_path"]
    os.environ["JULIA_PYTHONCALL_EXE"] = sys.executable
    os.environ["JULIA_CONDAPKG_BACKEND"] = "Null"
    os.environ["PYTHON_JULIACALL_STARTUP_FILE"] = "no"
    os.environ["PYTHON_JULIACALL_THREADS"] = "1"
    import juliacall

    jl = juliacall.newmodule("Pore2ProcessBridge")
    jl.seval("using Pore2ProcessDemo; import PythonCall")
    if str(jl.seval("Base.active_project()")) != str(ROOT / "julia/Project.toml"):
        raise RuntimeError("Julia loaded an unexpected project.")
    if str(jl.seval("pkgversion(PythonCall)")) != version("juliacall"):
        raise RuntimeError("PythonCall does not match JuliaCall. Rerun bootstrap.")
    return jl


def runtime_versions():
    _check_thread()
    jl = runtime()
    return {"python": platform.python_version(), "julia": str(jl.seval("VERSION")),
            **{name: version(name) for name in ("numpy", "panel", "pyvista", "vtk", "matplotlib", "juliacall")},
            "PythonCall": str(jl.seval("pkgversion(PythonCall)"))}


def run_case(inputs: NetworkInputs = NetworkInputs()) -> CaseResult:
    _check_thread()
    for name in ("nx", "ny", "nz"):
        value = getattr(inputs, name)
        if not isinstance(value, numbers.Integral) or isinstance(value, bool):
            raise ValueError(f"{name} must be an integer.")
    was_initialized = bool(runtime.cache_info().currsize)
    started = perf_counter()
    jl = runtime()
    initialized = perf_counter()
    try:
        s = jl.Pore2ProcessDemo.run_case(*asdict(inputs).values())
    except Exception as e:
        raise ValueError(f"Network input error: {_input_error(e)}") from e
    arrays = {}
    index_names = ("edges", "inlet_nodes", "outlet_nodes", "interior_nodes")
    array_names = ("coordinates_m", "edges", "pressures_Pa", "flows_m3_s", "radii_m", "lengths_m",
                   "conductances_m3_Pa_s", "net_outflows_m3_s", "inlet_nodes", "outlet_nodes", "interior_nodes")
    for name in array_names:
        values = np.array(getattr(s, name), dtype=np.int64 if name in index_names else np.float64, copy=True)
        if name in index_names:
            values -= 1  # The sole Julia -> Python node-index conversion.
        values.setflags(write=False)
        arrays[name] = values
    scalar_names = ("Q_in_m3_s", "Q_out_m3_s", "length_m", "area_m2", "max_interior_imbalance_m3_s",
                    "flow_mismatch_m3_s", "normalized_conservation_error", "linear_solve_s", "total_s")
    scalars = {name: float(getattr(s, name)) for name in scalar_names}
    k = None if s.permeability_m2 is None else float(s.permeability_m2)
    return CaseResult(inputs=inputs, run_id=uuid4().hex, timestamp_utc=datetime.now(timezone.utc).isoformat(),
                      **arrays, **scalars, permeability_m2=k,
                      runtime_initialization_s=0.0 if was_initialized else initialized-started,
                      bridge_call_s=perf_counter()-initialized)


def predict_column_flow(result: CaseResult, column: ColumnInputs = ColumnInputs()) -> float | None:
    _check_thread()
    try:
        value = runtime().Pore2ProcessDemo.predict_column_flow(result.permeability_m2, result.inputs.viscosity_Pa_s,
                                                              column.length_m, column.area_m2, column.delta_p_Pa)
    except Exception as e:
        raise ValueError(f"Column input error: {_input_error(e)}") from e
    return None if value is None else float(value)
