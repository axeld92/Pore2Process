"""Compare the live carrier with an independently relaxed, frozen salt plume."""
import argparse
import os
import shutil
import subprocess

from examples.adsorption.breakthrough_hypotheses.hydrodynamics import handoff
from .run import HERE, OUT


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case"); p.add_argument("label",choices=("early","hour"))
    p.add_argument("--threads",type=int,default=4)
    a=p.parse_args()
    parent=OUT/a.case
    folder=parent/("audit_"+a.label)
    folder.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(parent/"active.i32",folder/"active.i32")
    with (folder/"run.log").open("w") as log:
        subprocess.run([*handoff.julia(1),HERE/"frozen.jl",parent,a.label],cwd=folder,
                       env={**os.environ,"OMP_NUM_THREADS":str(a.threads),"OPENBLAS_NUM_THREADS":"1"},
                       stdout=log,stderr=subprocess.STDOUT,check=True)


if __name__=="__main__":
    main()
