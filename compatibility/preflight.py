"""Standard-library-only, bounded inventory of the tools used by these probes."""
from datetime import datetime, timezone
from pathlib import Path
import json
import os
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / ".local/compatibility"


def probe(executable, *arguments):
    path = shutil.which(str(executable))
    if path is None:
        return {"status": "absent", "requested": str(executable)}
    try:
        result = subprocess.run([path, *arguments], text=True, capture_output=True,
                                timeout=15, env={**os.environ, "DOTNET_CLI_HOME": str(CACHE / "dotnet-home"),
                                                "DOTNET_CLI_TELEMETRY_OPTOUT": "1",
                                                "DOTNET_GENERATE_ASPNET_CERTIFICATE": "false"})
        return {"path": path, "exit_code": result.returncode,
                "output": (result.stdout + result.stderr).strip()[:5000]}
    except subprocess.TimeoutExpired:
        return {"path": path, "status": "timed out after 15 s"}


def main():
    output = ROOT / "results/compatibility"
    output.mkdir(parents=True, exist_ok=True)
    records = {"timestamp_utc": datetime.now(timezone.utc).isoformat(),
               "project": str(ROOT), "os": platform.platform(), "machine": platform.machine(),
               "python": {"executable": sys.executable, "version": platform.python_version()},
               "tools": {name: probe(name, *arguments) for name, arguments in {
                   "g++-15": ["--version"], "gcc-15": ["--version"], "clang++": ["--version"],
                   "make": ["--version"], "cmake": ["--version"], "mpirun": ["--version"],
                   "git": ["--version"], "dotnet": ["--version"], "mono": ["--version"],
                   str(CACHE / "dotnet/dotnet"): ["--info"],
                   "pkg-config": ["--modversion", "tinyxml2", "zlib"],
               }.items()},
               "distributions": {},
               "native_dependencies": probe("otool", "-L", str(CACHE / "openlb_example_build/poiseuille2d"))}
    for name in ("distributions", "dotnet_distribution"):
        path = CACHE / (name + ".json")
        records["distributions"][name] = json.loads(path.read_text()) if path.exists() else None
    target = output / "environment.json"
    target.write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
