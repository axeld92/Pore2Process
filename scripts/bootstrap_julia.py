#!/usr/bin/env python3
"""Install/check exactly matching PythonCall in the single julia/ project."""
import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default=os.environ.get("PYTHON_JULIACALL_EXE") or shutil.which("julia"))
    args = parser.parse_args()
    if not args.julia:
        raise SystemExit("Julia not found. Install Julia via https://julialang.org/install/ then pass --julia /absolute/path/to/julia.")
    requested = version("juliacall")
    executable = subprocess.check_output([args.julia, "--startup-file=no", "-e", 'print(joinpath(Sys.BINDIR,"julia"))'], text=True).strip()
    project = ROOT / "julia"
    manifest = project / "Manifest.toml"
    installed = None
    if manifest.exists():
        data = tomllib.loads(manifest.read_text())
        installed = data.get("deps", {}).get("PythonCall", [{}])[0].get("version")
    setup = "using Pkg; "
    if installed != requested:
        setup += f'Pkg.add(PackageSpec(name="PythonCall",version=v"{requested}")); Pkg.compat("PythonCall","={requested}"); '
    setup += 'Pkg.instantiate(); Pkg.precompile(); using Pore2ProcessDemo; import PythonCall; '
    setup += f'@assert pkgversion(PythonCall) == v"{requested}"; '
    setup += f'@assert Base.active_project() == {json.dumps(str(project / "Project.toml"))}; println("Julia project verified: ",VERSION,"; PythonCall ",pkgversion(PythonCall))'
    env = os.environ.copy()
    env.update(JULIA_PYTHONCALL_EXE=sys.executable, JULIA_CONDAPKG_BACKEND="Null", JULIA_NUM_PRECOMPILE_TASKS="2")
    subprocess.run([executable, "--startup-file=no", f"--project={project}", "-e", setup], env=env, check=True)
    julia_version = subprocess.check_output([executable, "--startup-file=no", "-e", "print(VERSION)"], text=True).strip()
    config = {"executable": executable, "project": str(project), "julia_version": julia_version,
              "pythoncall_version": requested, "depot_path": env.get("JULIA_DEPOT_PATH")}
    target = ROOT / ".local/julia_runtime.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(config, indent=2) + "\n")
    print(f"Saved {target}; JuliaCall will use this executable and project, with no automatic Julia download.")


if __name__ == "__main__":
    main()
