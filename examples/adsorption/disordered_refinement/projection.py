"""Check the accelerated face projection and run the N96 handoff."""
import argparse
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.disordered_packing.mesh import FLOW_FILES, flow_checks, load_field, shared_faces
from examples.adsorption.disordered_packing.native import ORDERED
from examples.adsorption.dispersion_holdout.run import execute
from pore2process_demo.io import read_json, sha, write_json

from .native import GRID, HERE, PARENT, parent

INPUTS = (*FLOW_FILES,"coords.i32")
FILES = (*INPUTS,"balanced.f64")
PARITY_LIMIT = 1e-10


def project(directory,info,coords,faces):
    runtime = read_json(ORDERED/"frozen_prediction.json")["parity_command"][:4]
    target = np.array(info["shape"])*info["directional_flow_m3_s"]
    cmd = [*runtime,str(HERE/"balance.jl"),str(directory),str(len(coords)),*map(str,target)]
    print("Face projection start",directory,flush=True)
    timing = execute(cmd,directory/"balance.log",timeout=7200)
    result = np.fromfile(directory/"balanced.f64",dtype="<f8")
    checks = flow_checks(coords,info["shape"],info,faces,result)
    record = dict(info=info,checks=checks,elapsed=timing,command=cmd,
                  files={f:sha(directory/f) for f in FILES})
    print("Face projection complete",directory,checks,flush=True)
    return record,result


def parity(base):
    old = PARENT/"disordered_N64/faces"
    record = parent()
    for name in (*FILES,"summary.json"):
        if sha(old/name)!=record["artifacts_sha256"][f"disordered_N64/faces/{name}"]:
            raise ValueError("Pinned N64 projection evidence changed")
    directory = base/"projection_N64_parity"
    directory.mkdir()
    for name in INPUTS:
        shutil.copy2(old/name,directory/name)
    coords = np.fromfile(directory/"coords.i32",dtype="<i4").reshape(-1,3)
    faces = {name:np.fromfile(directory/name,dtype=dtype) for name,dtype in FLOW_FILES.items()}
    saved = read_json(old/"summary.json")
    result, q = project(directory,saved["info"],coords,faces)
    reference = np.fromfile(old/"balanced.f64",dtype="<f8")
    relative = float(np.max(abs(reference-q))/np.max(abs(reference)))
    if relative > PARITY_LIMIT:
        raise ValueError(f"Accelerated projection differs from the N64 reference by {relative:g}")
    result.update(max_face_difference_relative=relative,parity_limit=PARITY_LIMIT,
                  source_sha256=sha(HERE/"balance.jl"),reference_sha256=sha(old/"balanced.f64"))
    write_json(directory/"summary.json",result)
    print("N64 projection parity",relative,flush=True)


def handoff(base):
    verification = read_json(base/"projection_N64_parity/summary.json")
    if (verification["max_face_difference_relative"]>PARITY_LIMIT
            or verification["source_sha256"]!=sha(HERE/"balance.jl")):
        raise ValueError("Accelerated projection requires verified N64 parity")
    folder,coords,velocity,info = load_field(base,"disordered",GRID)
    directory = folder/"faces"
    directory.mkdir()
    faces = shared_faces(coords,velocity,info["shape"],info["raw"]["dx_m"])
    coords.astype("<i4").tofile(directory/"coords.i32")
    for name,values in faces.items():
        values.tofile(directory/name)
    result,_ = project(directory,info,coords,faces)
    write_json(directory/"summary.json",result)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base",type=Path)
    p.add_argument("--parity",action="store_true")
    args = p.parse_args()
    (parity if args.parity else handoff)(args.base.resolve())
