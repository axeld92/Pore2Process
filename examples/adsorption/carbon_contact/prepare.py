"""Run contact-instrumented trajectories on the frozen carrier field."""
from datetime import datetime, timezone
from uuid import uuid4

from examples.adsorption.dispersion_holdout.run import execute
from examples.adsorption.exchange_attribution.markers import MARKERS
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = ROOT/"examples/adsorption/carbon_contact"
OUTPUT = ROOT/"results/adsorption/carbon_contact"
ATTRIBUTION = ROOT/"results/adsorption/exchange_attribution/runs/20261112T095114Z-99e65489"
ATTRIBUTION_SHA = "a926d1471756ab8675685639e15b16c455f61aab5961576af96b7f6e9df803a8"


def specifications(base):
    if sha(ATTRIBUTION/"summary.json") != ATTRIBUTION_SHA:
        raise ValueError("Fixed-inventory sensitivity archive changed")
    attribution = validate_record(ATTRIBUTION/"summary.json")
    ref = validate_record(MARKERS/"summary.json")
    matching = read_json(ATTRIBUTION/"marker_plan.json")["cases"]
    sources = {f"reference_{k}":(ref["requests"][k], MARKERS, k, ref["artifacts_sha256"][k+".f64"])
               for k in ("primary","half_step","replicate")}
    sources.update({k:(item["command"],ATTRIBUTION,k,attribution["artifacts_sha256"][k+".f64"])
                    for k,item in matching.items()})
    results = {}
    for key,(cmd,archive,case,digest) in sources.items():
        original = archive/(case+".f64")
        command = [*cmd[:4],str(HERE/"track_contacts.jl"),cmd[5],str(base/(key+".f64")),*cmd[7:]]
        results[key] = dict(command=command,reference_pairs=str(original.relative_to(ROOT)),
                            reference_sha256=digest,diffusivity_m2_s=float(cmd[8]),dt_s=float(cmd[9]),count=int(cmd[11]))
    return results


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    cases = specifications(base)
    files = (HERE/"track_contacts.jl",HERE/"protocol.md",HERE/"prepare.py")
    write_json(base/"plan.json",dict(cases=cases,sources={str(p.relative_to(ROOT)):sha(p) for p in files}))
    print("Carbon contact marker cohort",base,flush=True)
    for name, spec in cases.items():
        elapsed = execute(spec["command"],base/(name+".log"),timeout=1800)
        write_json(base/(name+"_execution.json"),elapsed)
        print(name,elapsed,flush=True)
    print("Completed marker contact trajectories",base,flush=True)


if __name__=="__main__":
    main()
