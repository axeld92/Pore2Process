"""Run Julia without a browser and export a complete result bundle."""
import argparse
from dataclasses import fields
from pathlib import Path

from .bridge import NetworkInputs, run_case
from .export import write_exports
from .units import microlitres_per_minute


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("results/default"))
    for field in fields(NetworkInputs):
        parser.add_argument("--" + field.name.replace("_","-"),dest=field.name,
                            type=int if field.name in ("nx","ny","nz") else float, default=field.default)
    args = parser.parse_args()
    try:
        result = run_case(NetworkInputs(**{field.name:getattr(args,field.name) for field in fields(NetworkInputs)}))
        output = write_exports(result,args.output)
    except (ValueError, RuntimeError) as e:
        parser.exit(2,f"{e}\n")
    print(f"Julia run {result.run_id}: {len(result.pressures_Pa)} nodes / {len(result.edges)} throats")
    print(f"Inlet flow: {result.Q_in_m3_s:.15g} m^3/s ({microlitres_per_minute(result.Q_in_m3_s):.9g} µL/min)")
    print(f"Permeability: {result.permeability_m2} m^2; normalized conservation error: {result.normalized_conservation_error:.3g}")
    print(f"Exports: {output.resolve()}")


if __name__ == "__main__":
    main()
