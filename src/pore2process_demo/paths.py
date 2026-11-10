"""Project paths independent of example nesting and the current directory."""
import os
import json
from functools import lru_cache
from pathlib import Path

def project_root():
    configured = os.environ.get("PORE2PROCESS_ROOT")
    candidates = [Path(configured).expanduser().resolve()] if configured else Path(__file__).resolve().parents
    for candidate in candidates:
        if (candidate / "pyproject.toml").is_file() and (candidate / "julia/Project.toml").is_file():
            return candidate
    raise RuntimeError("Cannot find the Pore2Process checkout; set PORE2PROCESS_ROOT to its directory.")

ROOT = project_root()


@lru_cache(maxsize=4)
def _path_map(path):
    return json.loads(path.read_text())["paths"]


def project_path(path):
    """Resolve explicit project aliases; relative inputs are project-relative.

    External absolute paths and saved run directories retain their locations.
    """
    path = Path(path)
    if path.is_absolute():
        try:
            path = path.relative_to(ROOT)
        except ValueError:
            return path
    relative = path.as_posix()
    for name in ("reorganization.json", "resource_organization.json"):
        manifest = ROOT / "docs/history" / name
        if not manifest.is_file():
            continue
        for old, new in sorted(_path_map(manifest).items(), key=lambda item: -len(item[0])):
            if relative == old or relative.startswith(old + "/"):
                relative = new + relative[len(old):]
                break
    return ROOT / relative
