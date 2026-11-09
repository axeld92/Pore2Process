"""Build a byte-identical stock example in an external build directory."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / ".local/compatibility/openlb/olb-1.9r0"
BUILD = ROOT / ".local/compatibility/openlb_example_build"
OUTPUT = ROOT / "results/compatibility"


def main():
    # The shared external builder prepares a core archive without vendor writes.
    subprocess.run([sys.executable,str(ROOT / "hydraulic/build.py")],check=True,timeout=700)
    BUILD.mkdir(parents=True,exist_ok=True)
    OUTPUT.mkdir(parents=True,exist_ok=True)
    hashes={}
    for name in ("poiseuille2d.cpp","case.h"):
        source=SOURCE / "examples/laminar/poiseuille2d" / name
        target=BUILD / name
        if not target.exists() or target.read_bytes()!=source.read_bytes():
            shutil.copy2(source,target)
        hashes[name]=hashlib.sha256(source.read_bytes()).hexdigest()
    config=(ROOT / "hydraulic/openlb/config.mk").read_text()
    (BUILD / "config.mk").write_text(config)
    (BUILD / "Makefile").write_text("EXAMPLE = poiseuille2d\ninclude config.mk\ninclude $(OLB_ROOT)/default.single.mk\n")
    core=BUILD / "core"
    if not core.exists():
        core.symlink_to(ROOT / ".local/hydraulic/build/core",target_is_directory=True)
    command=["make","-j2","GIT_VERSION=1.9r0",f"OLB_ROOT={SOURCE}","onlysample"]
    started=perf_counter()
    with (OUTPUT / "openlb_build.log").open("w") as log:
        subprocess.run(command,cwd=BUILD,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    record={"version":"1.9r0","original_source":str(SOURCE),"external_build":str(BUILD),
            "build_command":command,"original_example_sha256":hashes,"config_mk":config,
            "build_s":perf_counter()-started,
            "compiler":subprocess.check_output(["/opt/homebrew/bin/g++-15","--version"],text=True).splitlines()[0]}
    (OUTPUT / "openlb_build.json").write_text(json.dumps(record,indent=2)+"\n")
    print(f"Built byte-identical poiseuille2d outside OpenLB: {BUILD}")


if __name__=="__main__":
    main()
