"""Write conservative Julia face fields from the rectangular OpenLB exports."""
import argparse
from math import gcd
from pathlib import Path
from time import sleep

import numpy as np

from examples.adsorption.disordered_packing import mesh
from examples.adsorption.dispersion_holdout.run import execute
from pore2process_demo.io import read_json, sha, write_json

from . import native


def project(base, flow, grid):
    folder, ijk, velocity, info = native.load_field(base, flow, grid)
    faces_dir = folder / "faces"
    summary = faces_dir / "summary.json"
    if faces_dir.exists():
        for _ in range(90):
            if summary.is_file():
                break
            sleep(2)
        else:
            raise ValueError(f"Previous face projection has no verified summary in {faces_dir}")
        _, _, _, saved = mesh.load_faces(folder)
        if saved["info"] != info:
            raise ValueError("Existing face projection belongs to a different field")
        return
    faces_dir.mkdir()
    shape = info["shape"]
    faces = mesh.shared_faces(ijk, velocity, shape, info["raw"]["dx_m"])
    for name, field in faces.items():
        field.tofile(faces_dir / name)
    ijk.astype("<i4").tofile(faces_dir / "coords.i32")
    _, _, probes = native.bypass.inputs()
    runtime = read_json(probes / "plan.json")["cases"][flow]["command"][:4]
    total = np.array(shape) * info["directional_flow_m3_s"]
    command = [*runtime, str(native.PACKING / "balance.jl"), str(faces_dir),
               str(len(ijk)), *map(str, total)]
    elapsed = execute(command, faces_dir / "balance.log", timeout=7200)
    q = np.fromfile(faces_dir / "balanced.f64", dtype="<f8")
    checks = mesh.flow_checks(ijk, shape, info, faces, q)
    write_json(faces_dir / "summary.json", dict(info=info, checks=checks, elapsed=elapsed,
        command=command, files={name: sha(faces_dir / name)
                                for name in (*mesh.FLOW_FILES, "coords.i32", "balanced.f64")}))
    print(f"{flow} N{grid} face projection", checks, flush=True)


def find_common_cut(minima):
    width = gcd(*(len(array) for array in minima))
    index = np.arange(width)
    valid = np.ones(width, dtype=bool)
    for array in minima:
        cut = index * (len(array) // width)
        minimum = array[(cut - 1) % len(array)]
        valid &= (minimum > 0) & np.isfinite(minimum)
    candidates = index[valid]
    if len(candidates) == 0:
        return None
    return int(candidates[0]), width


def compact(base, flows, grids):
    cases = {}
    for grid in grids:
        for flow in flows:
            folder = base / f"{flow}_N{grid}"
            ijk, faces, q, record = mesh.load_faces(folder)
            nx = record["info"]["shape"][0]
            cases[folder] = mesh.forward_cuts(ijk, faces, q, nx)
    match = find_common_cut(list(cases.values()))
    if match is None:
        details = {folder.name: int(np.count_nonzero(np.isfinite(values) & (values > 0)))
                   for folder, values in cases.items()}
        raise ValueError(f"No common entirely forward cut. Forward planes by case: {details}")
    index, width = match
    write_json(base / "cut.json", dict(index=index, width=width,
                                      folders=[f.name for f in cases]))
    print("Common physical forward cut", index, "/", width, flush=True)
    for folder, array in cases.items():
        nx = len(array)
        cut = index * (nx // width)
        ijk = np.memmap(folder / "faces/coords.i32", mode="r", dtype="<i4").reshape(-1, 3)
        offset = np.sort((ijk[:, 0] - cut) % nx)
        counts = np.searchsorted(offset, np.arange(nx // native.PERIODS, nx + 1, nx // native.PERIODS))
        if counts[-1] != len(ijk):
            raise ValueError("Incomplete rectangular compact periods")
        d = folder / "compact"
        summary = d / "summary.json"
        if summary.is_file():
            previous = read_json(summary)
            if (previous["cut_index"], previous["period_end_counts"], previous["face_summary_sha256"]) == (
                    cut, counts.tolist(), sha(folder / "faces/summary.json")):
                continue
        record = mesh.compact(folder, cut)
        record["period_end_counts"] = counts.tolist()
        record["files"] = {path.name: sha(path) for path in sorted(d.iterdir())
                           if path.suffix in (".i32", ".f64")}
        write_json(summary, record)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base", type=Path)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--project", type=int)
    mode.add_argument("--compact", type=int, nargs="+")
    p.add_argument("--flow", choices=native.FLOWS)
    args = p.parse_args()
    if args.project is not None:
        for flow in [args.flow] if args.flow else native.FLOWS:
            project(args.base.resolve(), flow, args.project)
    else:
        if args.flow:
            p.error("--flow applies only to --project")
        compact(args.base.resolve(), native.FLOWS, args.compact)


if __name__ == "__main__":
    main()
