"""Test the two literal Kd exponent sign conventions; never fit either dataset."""

from argparse import ArgumentParser
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from examples.biocatalysis.pentyl_propanoate.run import read_rows, write_rows
from examples.biocatalysis.pentyl_propanoate_equimolar import run as audit
from pore2process_demo.io import read_json,sha,write_json
from pore2process_demo.paths import ROOT

PROTOCOL = audit.HERE / "source_sign_protocol.json"
AUTHOR = {1:audit.SOURCE/"author_conversion_1to1.csv",2:audit.PUBLISHED.AUTHOR}
LABELS = ("literal_printed","opposite_exponent_sign")


def source_parameters(source, variant):
    values = np.array([source["rate_parameters"][k] for k in audit.PUBLISHED.PARAMETERS])
    index = audit.PUBLISHED.PARAMETERS.index("Kd_exp")
    if values[index] != -28091:
        raise ValueError("Unexpected Table S3 denaturation term.")
    if variant == "opposite_exponent_sign":
        values[index] = -values[index]
    elif variant != "literal_printed":
        raise ValueError(f"Unknown sign scenario: {variant}")
    return values


def author_curves():
    return {(feed,*audit.PUBLISHED.condition(row)):float(row[key])
            for (feed,file),key in zip(AUTHOR.items(),("author_predicted_conversion","acid_conversion"),strict=True)
            for row in read_rows(file)}


def predict_batch(jl, p, label, groups, cases, plan):
    rows=[]
    for assay,curve in groups.items():
        start=curve[0]
        feed,T,w=start["feed"],start["temperature_C"],start["loading_pct"]
        c=next(r for r in cases if r["feed_ratio"]==feed and r["temperature_C"]==T)
        # Actual charge-derived dose from the 2026 paper, not a fitted rate factor.
        C=(16.3 if feed==1 else 10.3)*.98807/.040*w/100
        fine_step=plan["predictions"]["batch_time_steps_min"][-1]
        for r in curve:
            t=r["time_min"]
            kwargs={"acid0":c["acid0_mol_L"],"alcohol0":c["alcohol0_mol_L"]}
            coarse,fine=[float(jl.LorenzoPackedBed.batch_prediction(p,T,t,C,**kwargs,max_step_min=dt))
                         for dt in plan["predictions"]["batch_time_steps_min"]]
            diff=abs(coarse-fine)
            if not (0<=fine<=1 and diff<=plan["predictions"]["batch_refinement_max_absolute_conversion"]):
                raise RuntimeError(f"Invalid batch replay for {label}, assay {assay}, t={t}, step={fine_step}")
            rows.append({"case":label,"assay":assay,"feed_ratio":feed,"temperature_C":T,
                         "loading_pct":w,"catalyst_g_L_nominal":C,"time_min":t,
                         "observed_conversion":r["printed"],"predicted_conversion":fine,
                         "batch_step_refinement_difference":diff})
    return rows


def run(output, previous_path):
    audit.check_prior()
    plan=read_json(PROTOCOL)
    main=read_json(audit.PLAN)
    source=read_json(audit.PUBLISHED.MODEL_INPUTS)
    jl=audit.load_runtime()
    cases=audit.conditions(jl,source,read_json(audit.EXTRA),read_json(audit.PROPERTIES))
    groups=audit.prior.read_batch(read_json(audit.prior.PLAN))
    batch,bed=[],[]
    for name in LABELS:
        print(f"Table S3 sign consistency: {name}",flush=True)
        p=jl.LorenzoPackedBed.PingPong(source_parameters(source,name))
        batch.extend(predict_batch(jl,p,name,groups,cases,plan))
        for c in cases:
            terms=jl.LorenzoPackedBed.rate_terms(p,c["temperature_C"],c["acid0_mol_L"],c["alcohol0_mol_L"],False)
            bed.append({"case":name,**audit.row_head(c),"Kd_numeric":float(terms.kd),
                        "Kx_from_rate_coefficients":float(terms.k1/terms.k2),
                        **audit.solve(jl,p,c,main["numerical_acceptance"])})
    if len(batch)!=450 or len(bed)!=80:
        raise ValueError("Incomplete two-sign source audit.")
    output.mkdir(parents=True,exist_ok=False)
    write_rows(output/"sign_bed_before_comparison.csv",bed)
    write_rows(output/"sign_batch_comparison.csv",batch)
    observed=audit.read_observations()
    authors=author_curves()
    if set(authors)!=set(observed):
        raise ValueError("Missing published model curve points.")
    cases_comp=[audit.add_observation(r,observed)|{"author_predicted_conversion":authors[audit.identity(r)]} for r in bed]
    write_rows(output/"sign_bed_comparison.csv",cases_comp)
    scores={}
    for name in LABELS:
        rows=[r for r in cases_comp if r["case"]==name]
        # Compare the same simulated conversions with the source's printed
        # calculated curves separately from actual reactor observations.
        author_scores=audit.score_groups([r|{"observed_conversion":r["author_predicted_conversion"]} for r in rows])
        batch_scores={}
        for f in (1,2):
            sel=[r for r in batch if r["case"]==name and r["feed_ratio"]==f and r["time_min"]>0]
            batch_scores[f"{f}to1"]=audit.scored(sel)
            hot=[r for r in sel if r["temperature_C"]>=70]
            batch_scores[f"{f}to1_70to80C"]=audit.scored(hot)
        scores[name]={"batch":batch_scores,"bed":audit.score_groups(rows),"author_curve":author_scores}
    # Prior literal model must be reproducible for the whole current matrix.
    previous={audit.identity(r):float(r["predicted_conversion"]) for r in read_rows(previous_path)
              if r["case"]==audit.REFERENCE_CASE}
    compare=max(abs(previous[audit.identity(r)]-r["predicted_conversion"]) for r in bed if r["case"]==LABELS[0])
    if set(previous)!=set(observed) or compare>1e-10:
        raise RuntimeError("Literal sign changes the earlier source replay.")
    author_only=[{"feed_ratio":k[0],"temperature_C":k[1],"flow_ml_min":k[2],
                  "predicted_conversion":val,"observed_conversion":observed[k]["observed_conversion"]} for k,val in authors.items()]
    paths=[PROTOCOL,audit.HERE/"source_sign_audit.py",audit.HERE/"run.py",audit.HERE/"model.jl",audit.PLAN,
           audit.PUBLISHED.MODEL_INPUTS,audit.PUBLISHED.MODEL,audit.EXTRA,audit.PROPERTIES,
           audit.prior.HERE/"run.py",audit.prior.PLAN,audit.prior.loading.film.HERE/"model.jl",
           *audit.prior.BATCH_FILES,*AUTHOR.values(),*audit.OBSERVATIONS.values()]
    summary={"scope":"Conditional sign consistency, not a new fit or a proven typo correction.",
             "julia":str(jl.seval("VERSION")),"batch_count":450,"bed_count":80,"catalyst_mass_kg":.006,
             "prior_literal_replay_absolute_conversion_difference":compare,
             "max_batch_step_refinement":max(r["batch_step_refinement_difference"] for r in batch),
             "max_bed_grid_change":max(r["grid_change_absolute_conversion"] for r in bed),
             "max_balance_error":max(r["normalized_balance_error"] for r in bed),
             "max_film_root_error_mol_L_min":max(r["film_balance_error_mol_L_min"] for r in bed),
             "author_curves_vs_observations":audit.score_groups(author_only),"scores":scores,
             "source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in paths},
             "previous_bed_comparison_sha256":sha(previous_path),
             "output_sha256":{p.name:sha(p) for p in sorted(output.iterdir())}}
    write_json(output/"sign_summary.json",summary)
    return summary


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--previous",type=Path,required=True,help="The current dual-feed bed_comparison.csv")
    parser.add_argument("--output",type=Path,default=audit.OUTPUT/("sign_"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")))
    args=parser.parse_args()
    record=run(args.output.resolve(),args.previous.resolve())
    for name,metrics in record["scores"].items():
        print(name,[(f,round(metrics["batch"][f]["rmse_percentage_points"],2),
                     round(metrics["bed"][f]["rmse_percentage_points"],2),
                     round(metrics["author_curve"][f]["rmse_percentage_points"],2)) for f in ("1to1","2to1")])
    print("Output:",args.output,flush=True)
