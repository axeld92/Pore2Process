"""Keep compact verification evidence from the reproducible full campaign."""
from argparse import ArgumentParser
from pathlib import Path
from shutil import copyfile

from examples.biocatalysis.pentyl_propanoate.run import read_rows,write_rows
from pore2process_demo.io import read_json,sha,write_json
from examples.biocatalysis.pentyl_propanoate_composition.run import HERE,LAWS

LEDGERS=("predictions_before_comparison.csv","bed_comparison.csv","stationary_clock_screen.csv")
PRIMARY_PLOT_CASES=("ratio_reference","reversible_shared","acid","acid_decay","protected_decay",
                    "residual_decay","protected_residual_decay","acid_both","protected_both")


def abridge_case(case):
    # Retain parameters, training identity, objective, stationarity and
    # per-curve metrics without duplicating dense Q/Jacobian diagnostics
    # or all full-run prediction rows.
    return {**{k:v for k,v in case.items() if k!="fit"},
        "fit":{law:{k:v for k,v in row.items() if k not in (
            "start_sse","singular_values","point_weighted_scores",
            "max_accepted_steps","max_rejected_steps",
            "validation_max_accepted_steps","validation_max_rejected_steps")}
            for law,row in case["fit"].items()}}


def retain(run,output):
    summary=read_json(run/"summary.json")
    b=summary["batch"]
    output.mkdir(parents=True,exist_ok=False)
    for name in LEDGERS:
        copyfile(run/name,output/name)
    data=read_rows(run/"batch_predictions.csv")
    main_names={"all_data",*b["paired_holdouts"].keys(),"hot_acid_challenge"}
    selected=[r for r in data if r["variant"]=="primary" and r["case"] in main_names
              and r["law"] in PRIMARY_PLOT_CASES]
    write_rows(output/"selected_batch_predictions.csv",selected)
    # Keep complete scores, equilibrium provenance and fit parameters
    # for the primary cases; sensitivity cases retain withheld scores
    # and fitted shape parameters. Full numerical outputs remain hashed.
    compact=dict(summary)
    compact["batch"]={
        "all_data":abridge_case(b["all_data"]),
        "paired_holdouts":{name:abridge_case(f) for name,f in b["paired_holdouts"].items()},
        "paired_scores":b["paired_scores"],"hot_acid_challenge":abridge_case(b["hot_acid_challenge"]),
        "active_domain_scores":b["active_domain_scores"],"qualification":b["qualification"],
        "critical_sensitivity_scores":{
            name:{case: {law:{"equal_assay_rmse_pp":v["equal_assay_rmse_pp"],
                              "scores_by_assay":v["scores_by_assay"],"parameters":v["parameters"],
                              "converged":v["converged"],"bounds":v["bounds"]}
                         for law,v in f["fit"].items()} for case,f in cases.items()}
            for name,cases in b["sensitivities"].items()},
        "exponent_sensitivity_scores":{
            name:{case:{law:{"equal_assay_rmse_pp":v["equal_assay_rmse_pp"],
                             "scores_by_assay":v["scores_by_assay"],"parameters":v["parameters"],
                             "bounds":v["bounds"],"converged":v["converged"]}
                        for law,v in f["fit"].items()} for case,f in cases.items()}
            for name,cases in b["exponent_sensitivities"].items()},
        "fixed_protection_paired_sensitivity":{k:{"scores":v["scores"],
            "per_fold_parameters":{case:{law:f["fit"][law]["parameters"] for law in
                ("protected_decay","protected_residual_decay")} for case,f in v["folds"].items()}}
            for k,v in b["fixed_protection_paired_sensitivity"].items()},
    }
    compact["retained_sha256"]={p.name:sha(p) for p in sorted(output.iterdir())}
    compact["full_summary_sha256"]=sha(run/"summary.json")
    write_json(output/"summary.json",compact)


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("--run",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    retain(a.run.resolve(),a.output.resolve())
