"""Bound early escape with a zero-concentration carbon surface."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.disordered_early_bypass.run import film_stages
from examples.adsorption.disordered_memory_column.run import transport
from pore2process_demo.io import read_json, sha, write_json


def escape(paths, total_ntu, repeats, allocation):
    """Propagate the steady clean-film flux through mixed eight-section repeats."""
    n = len(paths.flow)
    mass = np.zeros(n)
    mass[0] = 1.
    for donor, receiver, _, probability in film_stages(paths, total_ntu, 8*repeats, allocation):
        mass = np.bincount(receiver, weights=mass[donor]*probability, minlength=n)
    return float(mass.sum()**repeats)


def analyse(base):
    plan = base / "plan.json"
    if sha(plan) != read_json(base / "frozen.json")["plan_sha256"]:
        raise ValueError("Expected the frozen column plan")
    cases = {}
    for size, condition in read_json(plan)["cases"].items():
        p, access = condition["parameters"], condition["access"]
        paths = transport(size + "_primary", 4)
        area = 6*p["dry_mass_kg"]/(p["particle_density_kg_m3"]*p["particle_diameter_m"])
        total_ntu = p["film_m_s"]*area/p["flow_m3_s"]
        repeats = p["blocks"]
        np.testing.assert_allclose(escape(paths, 0., repeats, "wall"), 1., rtol=1e-12)
        available = access["direct_fraction"] + (1-access["direct_fraction"])*access["slow_ratio"]
        section_count = 8*repeats
        required = section_count * np.expm1(-np.log(.05)/section_count)
        results = {}
        for label, ntu in (("fixed_access", available*total_ntu), ("full_access", total_ntu)):
            ideal = (1 + ntu/section_count)**(-section_count)
            np.testing.assert_allclose(escape(paths, ntu, repeats, "flow"), ideal, rtol=1e-12)
            allocated = escape(paths, ntu, repeats, "wall")
            if allocated < ideal - 1e-12:
                raise ValueError("Contact allocation violates the clean-film bound")
            results[label] = dict(transfer_number=ntu, wall_allocated_escape=allocated,
                flow_proportional_escape=ideal, minimum_extra_conductance_fraction=max(0., required/ntu-1))
        cases[size] = results
    return dict(plan_sha256=sha(plan), target_escape=.05, cases=cases)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    args = parser.parse_args()
    record = analyse(args.base)
    write_json(args.base / "clean_limit.json", record)
    for size, result in record["cases"].items():
        print(size, result)
