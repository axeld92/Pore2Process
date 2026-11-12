"""Repeat a borderline time-step comparison on an additional marker ensemble."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.repeated_geometry.predict import observe
from pore2process_demo.io import read_json, sha, write_json

from .packing import PERIODS


def plan(directory,size):
    original = read_json(directory/"tracer_plan.json")
    count = original["marker_count"]
    commands = {}
    for label in ("primary","half_step"):
        case = size+"_"+label
        cmd = original["cases"][case]["command"].copy()
        cmd[6] = str(directory/(case+"_extra.f64"))
        cmd[12] = str(original["seed"]+count)
        commands[case] = cmd
    return dict(commands=commands,additional_count=count,
                original_plan_sha256=sha(directory/"tracer_plan.json"))


def combined_observations(directory,size):
    frozen = read_json(directory/(size+"_extra_plan.json"))
    if frozen != plan(directory,size):
        raise ValueError("Additional tracer ensemble plan changed")
    values = {}
    for key in frozen["commands"]:
        count = frozen["additional_count"]
        paths = [directory/(key+suffix+".f64") for suffix in ("","_extra")]
        data = np.concatenate([np.fromfile(path,dtype="<f8").reshape(count,2*PERIODS) for path in paths])
        values[key] = observe(data)
    return values


def run(directory,size):
    frozen = plan(directory,size)
    path = directory/(size+"_extra_plan.json")
    write_json(path,frozen)
    print("Additional tracer ensemble frozen",sha(path),flush=True)

    def run_case(item):
        key,command = item
        elapsed = execute(command,directory/(key+"_extra.log"),timeout=3600)
        write_json(directory/(key+"_extra_execution.json"),elapsed)
        print(key,elapsed,flush=True)

    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run_case,frozen["commands"].items()))
    write_json(directory/(size+"_combined_observations.json"),combined_observations(directory,size))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("size",choices=("d800","d450"))
    a = p.parse_args()
    run(a.directory.resolve(),a.size)
