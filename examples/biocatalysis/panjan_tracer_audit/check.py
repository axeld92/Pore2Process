"""Verify the retained evidence and prior-impulse identity without any solve."""
import argparse
from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT
from .audit import HERE,REF,PRIOR,FIGURE,prior_kernel,toml


def check(require_run=False):
    r=read_json(REF/"summary.json")
    for key,root in (("source_sha256",ROOT),("code_sha256",ROOT),("output_sha256",HERE)):
        for name,digest in r[key].items():
            if sha(root/name)!=digest:
                raise ValueError(f"Changed retained audit evidence: {name}")
    _,m,old=prior_kernel()
    if m!=r["archived_A_impulse"] or old["A"]["model_parameters"]!=r["frozen_delay_parameters"]:
        raise ValueError("Prior A identity changed")
    if r["reaction_enabled"] or r["new_native_run"]:
        raise ValueError("This is a dimensional/reduced audit only")
    if r["any_R2_reduced_case_pass"]!=any(v["all_gates"] for k,v in r["prediction_scores"].items() if k.startswith("R2__")):
        raise ValueError("Inconsistent R2 screening")
    base=ROOT/r["run_directory"]
    if require_run and not base.is_dir():raise FileNotFoundError(base)
    if base.is_dir():
        manifest=base/"prediction_manifest.json"; arrays=base/"predictions_before_comparison.npz"
        if sha(manifest)!=r["saved_prediction_manifest_sha256"] or sha(arrays)!=r["saved_prediction_sha256"]:
            raise ValueError("Saved forward prediction changed")
        if read_json(manifest)["prediction_sha256"]!=sha(arrays):
            raise ValueError("Forward-prediction manifest is inconsistent")
    pdf=ROOT/"literature/inbox/PanjanThesis.pdf"
    if pdf.is_file() and sha(pdf)!=toml(FIGURE)["pdf_sha256"]:
        raise ValueError("Private source PDF changed")
    return r


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--require-run",action="store_true")
    r=check(p.parse_args().require_run)
    print("Panjan audit hashes, prior TRT impulse and frozen A parameters agree; reaction:",r["reaction_enabled"])
