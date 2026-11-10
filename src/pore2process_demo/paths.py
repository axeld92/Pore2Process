"""Project paths independent of example nesting and the current directory."""
import os
from pathlib import Path

def project_root():
    configured = os.environ.get("PORE2PROCESS_ROOT")
    candidates = [Path(configured).expanduser().resolve()] if configured else Path(__file__).resolve().parents
    for candidate in candidates:
        if (candidate / "pyproject.toml").is_file() and (candidate / "julia/Project.toml").is_file():
            return candidate
    raise RuntimeError("Cannot find the Pore2Process checkout; set PORE2PROCESS_ROOT to its directory.")

ROOT = project_root()
