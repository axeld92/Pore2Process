"""Audit phenol retention from Figure 7 and the frozen N160 column inventories."""
import argparse
from pathlib import Path

import numpy as np

from examples.adsorption.alhamed.source import PDF, PDF_SHA, points
from examples.adsorption.alhamed.study import crossing
from examples.adsorption.disordered_repeat_join.run import JOINS, OUTPUT as JOIN_OUTPUT
from examples.adsorption.exchange_column.run import CASES, PIXELS, SOURCE_SHA, source_inputs
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/figure7_retention/summary.json"
REFERENCE = HERE/"reference_summary.json"
JOIN_RUN = JOIN_OUTPUT/"runs/20261116T160738Z-2328c4c8/summary.json"
JOIN_SHA = "5898c180e06eeebb2d00e32836bc91ccecf49b6d0679c95419476fa620cb0dae"


def figure_curve(markers):
    time = np.r_[0.,markers["time_s"]]
    fraction = np.r_[0.,markers["fraction"]]
    if not np.all(np.diff(time)>0) or not np.all((0<=fraction)&(fraction<=1)):
        raise ValueError("Expected ordered Figure 7 markers with concentrations between zero and feed")
    return time,fraction


def integrate(time, fraction, feed_rate, stop):
    """Trapezoidal integral through a query strictly within the observed window."""
    if not 0<=stop<=time[-1]:
        raise ValueError("Integration end must be inside the digitized window")
    nodes = np.r_[time[time<stop],stop]
    return float(feed_rate*np.trapezoid(1-np.interp(nodes,time,fraction),nodes))


def checked_model(run, basis):
    p,r = run["inputs"],run["result"]
    for key in ("dry_mass_kg","liquid_volume_m3","flow_m3_s","inlet_concentration_kg_m3",
                "capacity_kg_kg","affinity_m3_kg"):
        if p[key]!=basis[key]:
            raise ValueError(f"Frozen column input differs from the stated basis: {key}")
    time,out,liquid,solid,discharged = (
        np.array(r[k]) for k in ("time_s","outlet_concentration_kg_m3","liquid_inventory_kg",
                                 "adsorbed_inventory_kg","cumulative_outlet_kg"))
    feed_rate = basis["feed_rate_kg_s"]
    # The solver saves an interval mean at each right-hand time node.
    exported = np.r_[0.,np.cumsum(np.diff(time)*basis["flow_m3_s"]*out[1:])]
    np.testing.assert_allclose(exported,discharged,rtol=1e-10,atol=1e-14)
    np.testing.assert_allclose(liquid+solid,time*feed_rate-discharged,rtol=1e-10,atol=1e-14)
    return dict(time_s=time,outlet_fraction=out/basis["inlet_concentration_kg_m3"],
                retained_kg=liquid+solid,liquid_kg=liquid,adsorbed_kg=solid)


def model_snapshot(model, stop, feed_rate):
    time = model["time_s"]
    if not 0<stop<=time[-1]:
        raise ValueError("Model history does not cover the requested comparison")
    fraction = float(model["outlet_fraction"][np.searchsorted(time,stop,side="left")])
    return dict(**{key:float(np.interp(stop,time,model[key])) for key in ("retained_kg","liquid_kg","adsorbed_kg")},
                outlet_fraction=fraction,net_retention_rate_kg_s=feed_rate*(1-fraction))


def snapshot(time, fraction, models, stop, basis):
    flow = basis["feed_rate_kg_s"]
    total = integrate(time,fraction,flow,stop)
    outlet = float(np.interp(stop,time,fraction))
    solid = [max(0.,total-basis["full_void_liquid_capacity_kg"]),total]
    return dict(time_s=float(stop),figure7_retained_kg=total,figure7_outlet_fraction=outlet,
                figure7_net_retention_rate_kg_s=(1-outlet)*flow,
                figure7_solid_bound_kg=solid,figure7_loading_bound_kg_kg=[v/basis["dry_mass_kg"] for v in solid],
                model={key:model_snapshot(m,stop,flow) for key,m in models.items()})


def sensitivity(time, fraction, models, figure, feed_rate):
    tails = []
    gaps = {join:[] for join in JOINS}
    for sx in (-1,1):
        for sy in (-1,1):
            times = np.r_[0.,np.maximum(0.,time[1:]+sx*figure["reading_sensitivity_time_s"])]
            deficit = np.r_[1.,1-np.clip(fraction[1:]+sy*figure["reading_sensitivity_fraction"],0.,1.)]
            observed = float(feed_rate*np.trapezoid(deficit,times))
            tails.append(observed)
            for name,model in models.items():
                gaps[name].append(float(np.interp(times[-1],model["time_s"],model["retained_kg"])-observed))
    if np.any(np.diff(fraction)<0):
        raise ValueError("The monotone interpolation envelope requires nondecreasing marker values")
    widths = np.diff(time)
    return dict(reading_scenarios_retention_kg=[min(tails),max(tails)],
                reading_scenarios_model_excess_kg={key:[min(v),max(v)] for key,v in gaps.items()},
                monotone_interpolation_retention_kg=[
                    float(feed_rate*np.dot(widths,1-fraction[1:])),
                    float(feed_rate*np.dot(widths,1-fraction[:-1])),
                ])


def compare_case(markers, models, common_stop, figure, basis):
    time,fraction = figure_curve(markers)
    feed_rate = basis["feed_rate_kg_s"]
    mid = crossing(time,fraction,.5)
    if mid is None:
        raise ValueError("The experimental midpoint is outside the digitized window")
    common = snapshot(time,fraction,models,common_stop,basis)
    end = snapshot(time,fraction,models,time[-1],basis)
    midpoint = snapshot(time,fraction,models,mid,basis)
    model_stats = {}
    for key,model in models.items():
        excess = end["model"][key]["retained_kg"]-end["figure7_retained_kg"]
        midpoint_excess = midpoint["model"][key]["retained_kg"]-midpoint["figure7_retained_kg"]
        matched = crossing(model["time_s"],model["retained_kg"],end["figure7_retained_kg"])
        if matched is None:
            raise ValueError("Model does not bracket the measured end-window retention")
        model_stats[key] = dict(excess_at_last_marker_kg=excess,
            excess_built_after_observed_midpoint_fraction=(excess-midpoint_excess)/excess,
            model_at_fig7_end_retention=dict(time_s=matched,**model_snapshot(model,matched,feed_rate)))
    total = end["figure7_retained_kg"]
    deficit = basis["equilibrium_total_inventory_kg"]-total
    tail_rate = end["figure7_net_retention_rate_kg_s"]
    return dict(particle_diameter_m=markers["diameter_m"],point_count=len(time)-1,
                experimental_midpoint=midpoint,common_window_end=common,last_marker=end,
                model_diagnostics=model_stats,**sensitivity(time,fraction,models,figure,feed_rate),
                figure7_fraction_of_equilibrium_total=total/basis["equilibrium_total_inventory_kg"],
                unobserved_retention_to_reach_equilibrium_kg=deficit,
                minimum_monotone_tail_to_equilibrium_s=deficit/tail_rate)


def inputs():
    if sha(JOIN_RUN)!=JOIN_SHA or sha(PDF)!=PDF_SHA:
        raise ValueError("Frozen column run or published PDF changed")
    record = validate_record(JOIN_RUN)
    manifest = JOIN_RUN.parent/"frozen_predictions.json"
    if sha(manifest)!=record["frozen_prediction_sha256"] or record["figure7_fitted_parameters"]!=0:
        raise ValueError("Expected the frozen unfitted column predictions")
    frozen = read_json(manifest)
    _,figure,source,_ = source_inputs()
    if sha(PIXELS)!=figure["pixel_csv_sha256"]:
        raise ValueError("Figure 7 marker extraction changed")
    Q,C,M = (figure[k] for k in ("flow_m3_s","feed_kg_m3","derived_dry_mass_kg"))
    V = source["porosity"]*figure["length_m"]*figure["area_m2"]
    qmax,b = (source[k] for k in ("capacity_kg_kg","affinity_m3_kg"))
    qeq = qmax*(b*C)/(1+b*C)
    basis = dict(flow_m3_s=Q,inlet_concentration_kg_m3=C,feed_rate_kg_s=Q*C,dry_mass_kg=M,
        liquid_volume_m3=V,porosity=source["porosity"],bulk_density_kg_m3=source["bulk_density_kg_m3"],
        bed_length_m=figure["length_m"],bed_area_m2=figure["area_m2"],
        temperature_range_K=source["source_temperature_range_K"],capacity_kg_kg=qmax,affinity_m3_kg=b,
        equilibrium_loading_at_feed_kg_kg=qeq,equilibrium_adsorbed_inventory_kg=M*qeq,
        full_void_liquid_capacity_kg=V*C,equilibrium_total_inventory_kg=M*qeq+V*C)
    data = points()
    markers = {name:data[series] for name,(series,_) in CASES.items()}
    runs = {}
    for name,case in frozen["cases"].items():
        runs[name] = {}
        for join,item in case.items():
            curve = JOIN_RUN.parent/item["selected"]
            if sha(curve)!=frozen["curve_sha256"][curve.name]:
                raise ValueError("Frozen column history changed")
            runs[name][join] = checked_model(read_json(curve),basis)
    return figure,basis,markers,runs,record["frozen_prediction_sha256"]


def audit():
    figure,basis,markers,runs,digest = inputs()
    common = float(min(m["time_s"][-1] for m in markers.values()))
    return dict(schema="pore2process.figure7_retention.v1",
        source=dict(doi="10.1016/j.jhazmat.2009.05.002",figure=figure["location"],pdf_sha256=PDF_SHA,
                    pixel_sha256=figure["pixel_csv_sha256"],source_parameters_sha256=SOURCE_SHA,
                    reading_sensitivity_pixels=figure["reading_sensitivity_pixels"],pointwise_uncertainty=None,
                    per_fraction_weighed_carbon_mass_kg=None),
        model=dict(summary=JOIN_RUN.relative_to(ROOT).as_posix(),summary_sha256=JOIN_SHA,
                   frozen_prediction_sha256=digest,new_fitted_parameters=0),
        basis=basis,common_window_end_s=common,
        cases={name:compare_case(m,runs[name],common,figure,basis) for name,m in markers.items()})


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.check:
        if result!=read_json(REFERENCE):
            raise ValueError("Retention audit differs from the tracked reference")
        print("Retention audit matches the tracked reference")
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        write_json(args.output,result)
        print(args.output)
