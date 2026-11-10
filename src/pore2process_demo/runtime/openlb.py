"""Build the external case/core without writing into the OpenLB distribution."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import tarfile
from time import perf_counter

from pore2process_demo.paths import ROOT
OPENLB = ROOT / ".local/compatibility/openlb/olb-1.9r0"
BUILD = ROOT / ".local/hydraulic/build"
OUTPUT = ROOT / "results/hydraulic"


def verify_vendor():
    """Compare original source files with the pinned archive (no write operations)."""
    archive = ROOT / ".local/compatibility/downloads/olb-1.9r0.tgz"
    expected = "e4d2421b50643482036917c4416145f7d6b495492f3afac64e048212361c8b3e"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
        raise RuntimeError("OpenLB distribution archive hash changed.")
    checked = 0
    with tarfile.open(archive) as package:
        for entry in package:
            if not entry.isfile():
                continue
            relative = Path(entry.name).relative_to("olb-1.9r0")
            # config.mk was configured by the earlier compatibility probe.
            # It is never loaded or changed by this external build.
            if str(relative) == "config.mk":
                continue
            target = OPENLB / relative
            if target.read_bytes() != package.extractfile(entry).read():
                raise RuntimeError(f"Original OpenLB file differs from archive: {relative}")
            checked += 1
    return {"archive_sha256": expected, "original_files_checked": checked,
            "preexisting_vendor_config_ignored": True}


def build_case(source_directory=ROOT / "examples/hydraulics/capillary_benchmark/openlb", build_directory=BUILD, output_directory=OUTPUT):
    """Build an application folder with the same external, stock OpenLB core."""
    source_directory = Path(source_directory).resolve()
    BUILD, OUTPUT = Path(build_directory).resolve(), Path(output_directory).resolve()
    for path in (source_directory, BUILD, OUTPUT):
        if path == OPENLB.resolve() or OPENLB.resolve() in path.parents:
            raise ValueError("Application code, builds and results must stay outside original OpenLB.")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    vendor = verify_vendor()
    BUILD.mkdir(parents=True, exist_ok=True)
    core = BUILD / "core"
    core.mkdir(exist_ok=True)
    hashes = {}
    for source in sorted(source_directory.iterdir()):
        if not source.is_file():
            continue
        hashes[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
        target = BUILD / source.name
        if not target.exists() or target.read_bytes() != source.read_bytes():
            shutil.copy2(source, target)
    # These are the five translation units listed by the upstream core target.
    # Only object/archive locations differ; the original files are compiled directly.
    sources = ("communication/mpiManager.cpp", "communication/ompManager.cpp", "core/olbInit.cpp",
               "core/expr.cpp", "io/ostreamManager.cpp")
    compiler = "/opt/homebrew/bin/g++-15"
    flags = ["-O2", "-Wall", "-std=c++20", "-pthread", "-fPIC", "-DPLATFORM_CPU_SISD",
             '-DOLB_VERSION="1.9r0"', f"-I{OPENLB / 'src'}", "-I/opt/homebrew/include"]
    configuration={"flags":flags,"compiler":subprocess.check_output([compiler,"--version"],text=True),
                   "external_config":(BUILD / "config.mk").read_text(),
                   "vendor_archive":vendor["archive_sha256"]}
    stamp=BUILD / "configuration.json"
    changed=not stamp.exists() or json.loads(stamp.read_text())!=configuration
    started = perf_counter()
    commands = []
    with (OUTPUT / "build.log").open("w") as log:
        objects = []
        for relative in sources:
            source = OPENLB / "src" / relative
            target = core / (relative.replace("/", "_") + ".o")
            objects.append(target)
            if changed or not target.exists():
                command = [compiler, *flags, "-c", str(source), "-o", str(target)]
                commands.append(command)
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
        library = core / "libolbcore.a"
        if not library.exists() or any(p.stat().st_mtime > library.stat().st_mtime for p in objects):
            command = ["ar", "rcs", str(library), *map(str, objects)]
            commands.append(command)
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=30)
        if changed:
            (BUILD / "capillary3d.o").unlink(missing_ok=True)
        # The upstream onlysample target avoids the in-tree core/dependency targets.
        command = ["make", "-j2", "GIT_VERSION=1.9r0", f"OLB_ROOT={OPENLB}", f"PROJECT_ROOT={ROOT}", "onlysample"]
        commands.append(command)
        result = subprocess.run(command, cwd=BUILD, stdout=log, stderr=subprocess.STDOUT, timeout=600)
    if result.returncode:
        print((OUTPUT / "build.log").read_text()[-7000:])
        raise SystemExit(f"Build failed; see {OUTPUT / 'build.log'}")
    stamp.write_text(json.dumps(configuration,indent=2)+"\n")
    metadata = {"commands": commands, "source_sha256": hashes, "build_s": perf_counter()-started,
                "external_config": (BUILD / "config.mk").read_text(), "vendor_integrity": vendor,
                "core_sources": list(sources), "core_flags": flags,
                "compiler": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]}
    verify_vendor()
    (OUTPUT / "build.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(f"Built external capillary3d in {metadata['build_s']:.2f} s; {vendor['original_files_checked']} original files verified.")
    return metadata


def main():
    build_case()


if __name__ == "__main__":
    main()
