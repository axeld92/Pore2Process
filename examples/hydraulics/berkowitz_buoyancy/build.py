"""Build the external salt case with the same stock core as the saved carrier."""
from pathlib import Path
import subprocess

from examples.adsorption.breakthrough_hypotheses.hydrodynamics import native
from pore2process_demo.io import sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.runtime.openlb import OPENLB, verify_vendor

HERE = Path(__file__).resolve().parent
OUT = ROOT / ".local/berkowitz_buoyancy"
BINARY = OUT / "salt3d"


def build():
    stock = verify_vendor()
    core = native.BUILD / "build/core/libolbcore.a"
    if not core.is_file():
        native.build()
    src = HERE / "openlb/salt3d.cpp"
    OUT.mkdir(parents=True, exist_ok=True)
    compiler = "/opt/homebrew/bin/g++-15"
    command = [compiler, "-O3", "-march=native", "-std=c++20", "-pthread", "-fPIC",
               "-DPLATFORM_CPU_SISD", "-DPARALLEL_MODE_OMP", "-fopenmp", '-DOLB_VERSION="1.9r0"',
               f"-I{OPENLB/'src'}", "-I/opt/homebrew/include", str(src), str(core),
               "-L/opt/homebrew/lib", "-ltinyxml2", "-lz", "-o", str(BINARY)]
    with (OUT/"build.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    write_json(OUT/"build.json", dict(source_sha256=sha(src), binary_sha256=sha(BINARY),
               core_sha256=sha(core), vendor=stock, command=command))
    verify_vendor()


if __name__ == "__main__":
    build()
