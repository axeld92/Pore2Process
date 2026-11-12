"""Run first-hit clocks at the same two Péclet numbers and time steps."""
from datetime import datetime, timezone
from uuid import uuid4

from examples.adsorption.contact_column.inputs import CONTACT, CONTACT_SHA
from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = ROOT/"examples/adsorption/wall_uptake"
OUTPUT = ROOT/"results/adsorption/wall_uptake"


def requests(base):
    path = CONTACT/"summary.json"
    if sha(path)!=CONTACT_SHA:
        raise ValueError("Qualified contact cohort changed")
    archived = validate_record(path)
    plans = read_json(CONTACT/"plan.json")["cases"]
    cases = {}
    for name in ("d800_primary","d800_half_step","d450_primary","d450_half_step"):
        cmd = plans[name]["command"]
        cases[name] = dict(command=[*cmd[:4],str(HERE/"track_uptake.jl"),cmd[5],str(base/(name+".f64")),*cmd[7:]],
                           dt_s=plans[name]["dt_s"],diffusivity_m2_s=plans[name]["diffusivity_m2_s"],
                           reference_sha256=archived["artifacts_sha256"][name+".f64"])
    return cases


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    cases = requests(base)
    files = (HERE/"prepare.py",HERE/"track_uptake.jl",HERE/"protocol.md")
    write_json(base/"plan.json",dict(cases=cases,sources={str(f.relative_to(ROOT)):sha(f) for f in files}))
    print("Wall uptake cohort",base,flush=True)
    for name,item in cases.items():
        elapsed = execute(item["command"],base/(name+".log"),timeout=1800)
        write_json(base/(name+"_execution.json"),elapsed)
        print(name,elapsed,flush=True)
    print("Complete wall uptake marker tracking",base,flush=True)


if __name__=="__main__":
    main()
