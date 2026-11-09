#!/usr/bin/env python3
"""Bounded, standard-library-only inventory. No package installation or env dump."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ("numpy", "panel", "pyvista", "vtk", "matplotlib", "juliacall", "pytest")


def command(args):
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=20)
        return {"returncode": p.returncode, "output": (p.stdout + p.stderr).strip()[:2000]}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {"error": str(e)}


def inventory(https=False):
    report = {"project": str(ROOT), "os": platform.platform(), "architecture": platform.machine(),
              "python": {"executable": sys.executable, "version": platform.python_version()},
              "tools": {}, "packages": {}, "files": {}}
    for name in ("julia", "juliaup", "uv", "git", "codex"):
        path = shutil.which(name)
        if name == "uv" and path is None and (ROOT / ".tools/bin/uv").exists():
            path = str(ROOT / ".tools/bin/uv")
        report["tools"][name] = {"path": path, **(command([path, "--version"]) if path else {})}
    for name in PACKAGES:
        try:
            version = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            report["packages"][name] = {"status": "absent"}
            continue
        if name == "juliacall":
            # Never trigger JuliaPkg's automatic downloads just for an inventory.
            config = ROOT / ".local/julia_runtime.json"
            if not config.exists():
                report["packages"][name] = {"version": version, "status": "installed; import deferred until bootstrap"}
                continue
            probe = "from pore2process_demo.bridge import runtime; runtime(); print('import passed')"
        else:
            probe = f"import {name}; print('import passed')"
        result = command([sys.executable, "-c", probe])
        report["packages"][name] = {"version": version, "status": "import passed" if result.get("returncode") == 0 else "import failed", **result}
    for name in (".venv", "uv.lock", "julia/Project.toml", "julia/Manifest.toml", ".local/julia_runtime.json"):
        report["files"][name] = (ROOT / name).exists()
    report["git"] = command(["git", "status", "--short", "--branch"])
    if https:
        report["https"] = {}
        for url in ("https://pypi.org/pypi/panel/json", "https://pkg.julialang.org", "https://julialang-s3.julialang.org/bin/versions.json"):
            try:
                with urllib.request.urlopen(url, timeout=8) as response:
                    report["https"][url] = {"status": response.status}
            except Exception as e:
                report["https"][url] = {"error": str(e)}
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--https", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(inventory(args.https), indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)
