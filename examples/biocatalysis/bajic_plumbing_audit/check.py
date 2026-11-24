"""Verify the tube moment and/source reference, without claiming a new run."""
import argparse
from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT
from .audit import HERE,REF,BAJIC,STRNISA,PANJAN,toml


def check(require_run=False):
    s=read_json(REF/"summary.json")
    for section,root in (("source_sha256",ROOT),("code_sha256",ROOT),("output_sha256",HERE)):
        for key,digest in s[section].items():
            if sha(root/key)!=digest:raise ValueError(f"Changed retained evidence: {key}")
    if s["reaction_enabled"] or s["new_porous_OpenLB_run"] or not s["actual_new_Julia_tube_run"]:
        raise ValueError("Wrong numerical scope")
    if s["any_same_id_tube_necessary_variance_pass"]!=any(v["necessary_variance_screen_pass"] for r in s["residuals"] for v in r["diffusivity_cases"].values()):
        raise ValueError("Inconsistent necessary screen")
    if sha(REF/"tube_moments.csv")!=s["tube_run"]["table_sha256"]:raise ValueError("Julia table does not match")
    b=ROOT/s["run_directory"]
    if require_run and not b.is_dir():raise FileNotFoundError(b)
    if b.is_dir():
        if sha(b/"inputs.json")!=s["input_manifest_sha256"] or sha(b/"tube_moments.csv")!=s["tube_run"]["table_sha256"]:
            raise ValueError("Saved Julia run changed")
    for name,source in (("Bajic2017Development.pdf",BAJIC),("Strnisa2018Characterization.pdf",STRNISA),("PanjanThesis.pdf",PANJAN)):
        path=ROOT/"literature/inbox"/name
        if path.is_file() and sha(path)!=toml(source)["pdf_sha256"]:raise ValueError(f"Source PDF changed: {name}")
    return s


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--require-run",action="store_true")
    s=check(p.parse_args().require_run)
    print("Bajić plumbing evidence and actual Julia tube moment hashes agree; reaction:",s["reaction_enabled"])
