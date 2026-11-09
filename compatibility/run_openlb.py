"""Run unmodified OpenLB analytical channel benchmarks, with bounded acceptance."""
from datetime import datetime, timezone
from pathlib import Path
import json
import math
import re
import subprocess
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / ".local/compatibility/openlb_example_build/poiseuille2d"
OUTPUT = ROOT / "results/compatibility/openlb"


def main():
    cases = []
    # Binary-exact mesh spacings avoid an observed stock-example outlet selection
    # failure at N=40. Keep that exploratory failure documented; check materials.
    for resolution in (32, 64):
        directory = OUTPUT / f"n{resolution}"
        directory.mkdir(parents=True, exist_ok=True)
        arguments = {
            "RESOLUTION": resolution, "REYNOLDS": 10,
            "DOMAIN_EXTENT": "[2 1]", "PHYS_CHAR_VELOCITY": 1,
            "PHYS_CHAR_DENSITY": 1, "PHYS_STAT_ITER_T": 1.5,
            "OVERLAP": 3, "TUNER_PARAM": 0.5,
            "LATTICE_RELAXATION_TIME": 0.8, "MAX_PHYS_T": 30,
            "CONVERGENCE_PRECISION": 1e-9, "CONV_ITER": 0.25,
            # In this upstream example 0 actually enters writePNG(); 1 collects
            # EOC data without invoking gnuplot. Neither enables VTK graphics.
            "VTK_ENABLED": 0, "GNUPLOT_ENABLED": 1,
            "COMPUTE_ERROR": 1, "FLOW_TYPE": 1, "BOUNDARY_TYPE": 2,
        }
        command = [str(BINARY)] + [v for key, value in arguments.items() for v in ("--" + key, str(value))]
        started = perf_counter()
        with (directory / "stdout.log").open("w") as log:
            print(f"Running OpenLB channel benchmark N={resolution}...", flush=True)
            subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=300)
        elapsed = perf_counter() - started
        output = (directory / "stdout.log").read_text()
        errors = {}
        for field, norm, kind, value in re.findall(
                r"(velocity|pressure|strainRate)-(L1|L2|Linf)-error\((abs|rel)\)=([^\s;]+)", output):
            errors[f"{field}_{norm}_{kind}"] = float(value)
        iterations = re.findall(r"Converged after (\d+) iterations", output)
        materials = {int(material): int(count) for material, count in
                     re.findall(r"materialNumber=(\d+); count=(\d+)", output)}
        if len(errors) != 18 or not all(math.isfinite(value) for value in errors.values()):
            raise RuntimeError(f"Missing or non-finite analytical errors for N={resolution}; inspect stdout.log.")
        if any(materials.get(material) != resolution - 1 for material in (3, 4)):
            raise RuntimeError(f"N={resolution} inlet/outlet material counts incorrect: {materials}.")
        if not iterations or errors["velocity_L2_rel"] >= 0.01 or errors["pressure_L2_rel"] >= 0.05:
            raise RuntimeError(f"N={resolution} failed convergence/velocity/pressure acceptance; inspect stdout.log.")
        cases.append({"resolution": resolution, "arguments": arguments, "command": command,
                      "material_counts": materials,
                      "converged_iterations": int(iterations[-1]), "errors": errors,
                      "wall_s": elapsed, "log": str(directory / "stdout.log")})
        print(f"N={resolution}: converged; velocity relative L2 = {errors['velocity_L2_rel']:.6g}; {elapsed:.2f} s", flush=True)
    ratio = cases[0]["errors"]["velocity_L2_rel"] / cases[1]["errors"]["velocity_L2_rel"]
    pressure_ratio = cases[0]["errors"]["pressure_L2_rel"] / cases[1]["errors"]["pressure_L2_rel"]
    if ratio <= 1 or pressure_ratio <= 1:
        raise RuntimeError("Velocity/pressure error did not decrease under refinement.")
    result = {"tool": "OpenLB", "version": "1.9r0", "case": "official poiseuille2d",
              "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "binary": str(BINARY), "cases": cases, "velocity_L2_refinement_ratio": ratio,
              "pressure_L2_refinement_ratio": pressure_ratio,
              "scope": "2D channel compatibility/analytical benchmark, not porous-material permeability or 3D network validation"}
    (OUTPUT / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Passed: {OUTPUT / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
