"""Cold-process repeat of primary fits and frozen bed predictions."""
from argparse import ArgumentParser
from pathlib import Path
import numpy as np

from examples.biocatalysis.pentyl_propanoate_composition import run as r
from examples.biocatalysis.pentyl_propanoate_composition.check import check,REFERENCE
from pore2process_demo.io import read_json,write_json
from examples.biocatalysis.pentyl_propanoate.run import read_rows


def replay():
    check()
    expected=read_json(REFERENCE/"summary.json")
    plan=r.read_json(r.PLAN)
    groups=r.read_batch(plan)
    eq=r.equilibrium_rows()
    jl=r.load_runtime()
    model=jl.CompositionKinetics
    source=r.read_json(r.ACID_SOURCE)
    opts={"acid_exponent":source["inhibition_exponent"],
          "analogue_KI":source["recalculated_inhibition_concentration_mol_per_kg_liquid_substrates"]}
    ids=sorted(groups)
    cases={"all_data":(ids,ids)} | {name:([a for a in ids if a not in test],test)
        for name,test in plan["holdouts"]["paired_folds"].items()} | {
        "hot_acid_challenge":([a for a in ids if a<15],[15,16])}
    parameter_change=0.;score_change=0.;predictions=[]
    for name,(train,test) in cases.items():
        fit,rows=r.fit_case(model,groups,eq,name,"primary",train,test,**opts)
        ref=(expected["batch"]["all_data"] if name=="all_data" else
             expected["batch"]["hot_acid_challenge"] if name=="hot_acid_challenge" else
             expected["batch"]["paired_holdouts"][name])
        if fit["equilibrium_training"]!=ref["equilibrium_training"]:
            raise ValueError("Training equilibrium fit did not replay exactly.")
        for law in r.LAWS:
            got=fit["fit"][law];want=ref["fit"][law]
            p=np.array(list(got["parameters"].values()))
            q=np.array(list(want["parameters"].values()))
            parameter_change=max(parameter_change,float(np.max(abs(p-q)/np.maximum(1,abs(q)))))
            score_change=max(score_change,abs(got["equal_assay_rmse_pp"]-want["equal_assay_rmse_pp"]))
        predictions.extend(rows)
    key=lambda x:(x["case"],x["law"],int(x["assay"]),float(x["time_min"]))
    got={key(x):x for x in predictions}
    ref=read_rows(REFERENCE/"selected_batch_predictions.csv")
    progress_change=max(abs(got[key(x)]["predicted_conversion"]-float(x["predicted_conversion"])) for x in ref)
    activity_change=max(abs(got[key(x)]["relative_catalyst_activity"]-float(x["relative_catalyst_activity"])) for x in ref)
    bed,clock,_=r.bed_campaign(jl,expected["batch"],plan)
    bkey=lambda x:(x["case"],int(x["temperature_C"]),float(x["flow_ml_min"]))
    br={bkey(x):x for x in read_rows(REFERENCE/"predictions_before_comparison.csv")}
    bed_change=max(abs(x["predicted_conversion"]-float(br[bkey(x)]["predicted_conversion"])) for x in bed)
    ckey=lambda x:(x["law"],int(x["temperature_C"]),float(x["flow_ml_min"]),float(x["catalyst_age_min"]))
    cr={ckey(x):x for x in read_rows(REFERENCE/"stationary_clock_screen.csv")}
    clock_change=max(abs(x[k]-float(cr[ckey(x)][k])) for x in clock
        for k in ("conversion_uniform_inlet_activity","conversion_uniform_equilibrium_activity"))
    result={"primary_cases_repeated":len(cases),"models_per_case":len(r.LAWS),
        "selected_batch_points_compared":len(ref),"bed_points_compared":len(bed),"clock_points_compared":len(clock),
        "max_relative_parameter_change":parameter_change,"max_equal_curve_score_change_pp":score_change,
        "max_batch_conversion_change":progress_change,"max_batch_activity_change":activity_change,
        "max_bed_conversion_change":bed_change,"max_frozen_clock_conversion_change":clock_change}
    if max(parameter_change,score_change,progress_change,activity_change,bed_change,clock_change)>1e-9:
        raise RuntimeError(f"Primary reproducibility check failed: {result}")
    return result


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=replay()
    print(result)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        write_json(args.output.resolve(),result)
