"""Check retained evidence/holdout identity without executing a simulation."""
import argparse

from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT
from .reporting import REF,HERE,mean_delay_bound
from .compare import numerical_ok,passes


def check(require_run=False):
    r=read_json(REF/"summary.json")
    for section,parent in (("source_sha256",ROOT),("code_sha256",ROOT),("output_sha256",HERE)):
        for f,digest in r[section].items():
            if sha(parent/f)!=digest:
                raise ValueError(f"Changed retained evidence: {parent/f}")
    a=r["A"]
    c=r["C_manifest"]
    p=r["plan"]["protocol"]
    selected=a["selected"]["tracer"]
    held=c["carrier"]["tracer"]
    if r["reaction_enabled"] or a["reaction_enabled"] or c["reaction_enabled"]:
        raise ValueError("Reaction is outside this qualification")
    if (selected["reactor"],held["reactor"],a["model_parameters"]["calibration_reactor"])!=("A","C","A"):
        raise ValueError("Invalid calibration boundary")
    for key in ("depth_cells","dt_s","D_m2_s","factor","seed","particle_count"):
        if selected[key]!=held[key]:
            raise ValueError(f"Long tracer changed a frozen setting: {key}")
    if a["numerical_qualification"] != bool(a["grid_pass"] and a["time_step_pass"] and numerical_ok(a["selected"],p)):
        raise ValueError("Short numerical gate is inconsistent")
    if c["numerical_ok"]!=numerical_ok(c["carrier"],p):
        raise ValueError("Long numerical gate is inconsistent")
    for scores in (a["A_scores"],r["C_scores"]):
        for name,score in scores.items():
            gates=passes(score,p)
            if score["gates"]!=gates or score["all_gates"]!=all(gates.values()):
                raise ValueError(f"Measurement gate inconsistent: {name}")
    if a["independent_A_pass"]!=bool(a["numerical_qualification"] and a["A_scores"]["measured_blank_bg0"]["all_gates"]):
        raise ValueError("Short experimental gate is inconsistent")
    if r["independent_C_pass"]!=bool(c["numerical_ok"] and r["C_scores"]["measured_blank_bg0"]["all_gates"]):
        raise ValueError("Long experimental gate is inconsistent")
    if r["reaction_gate_pass"]!=bool(a["independent_A_pass"] and r["independent_C_pass"]):
        raise ValueError("Reaction gate is inconsistent")
    if (r["compact_record_sha256"]["C_prediction_before_comparison.npz"]!=c["predictions_sha256"] or
        r["compact_record_sha256"]["frozen_A.json"]!=c["frozen_A_sha256"]):
        raise ValueError("Retained holdout manifest is inconsistent")
    if (mean_delay_bound(a["A_scores"],r["C_scores"],r["plan"],a["solid_volume_ratio_C"])!=
        r["positive_delay_mean_screen"]):
        raise ValueError("Positive-delay mean feasibility calculation changed")
    base=ROOT/r["run_directory"]
    if require_run and not base.is_dir():
        raise FileNotFoundError(base)
    if base.is_dir():
        for f,digest in r["compact_record_sha256"].items():
            if sha(base/f)!=digest:
                raise ValueError(f"Changed saved run: {f}")
        if sha(base/"frozen_A.json")!=r["C_manifest"]["frozen_A_sha256"]:
            raise ValueError("Long prediction did not use frozen A")
    return r


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-run",action="store_true")
    args=parser.parse_args()
    r=check(args.require_run)
    print("Retained Strnisa tracer evidence matches its source/code/output hashes;",
          "reaction gate:",r["reaction_gate_pass"])
