"""Couple the frozen N160 tracer transport to the batch-qualified particles."""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from uuid import uuid4

import numpy as np

from examples.adsorption.alhamed.study import compare
from examples.adsorption.contact_column.run import crossings
from examples.adsorption.disordered_travel_memory.run import PROBES, inputs as tracer_inputs
from examples.adsorption.exchange_column.run import (
    BATCH_SHA, CASES, JULIA, PARTICLE, PIXELS, SOURCE, source_inputs, validate_solution,
)
from examples.adsorption.film_particle.bridge import call
from examples.adsorption.film_particle.metrics import column_passed
from examples.adsorption.shared_contact.run import full_refinement, solver as shared_solver
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/"results/adsorption/disordered_memory_column"
MEMORY = ROOT/"results/adsorption/disordered_travel_memory/runs/20261116T153803Z-b1f18981"
MEMORY_SHA = "5050b9b8a58682cf2902c1a3208d3cad3d2acff693d7f423fb6ac2619abfc734"
LEVELS = ((24,8.,.02),(48,4.,.01),(96,2.,.005),(128,1.,.0025))
CONTROLS = dict(flow_area=("primary",4,False),one_age_group=("primary",1,True),
                half_probe_step=("half_step",4,True))


@lru_cache(maxsize=1)
def julia_models():
    shared_solver()
    jl = runtime()
    for file in ("model.jl","handoff.jl"):
        jl.include(str(HERE/file))
    return jl.ColumnHandoff.load_transport, jl.MemoryColumn.solve_column


@lru_cache(maxsize=8)
def transport(name, groups):
    load, _ = julia_models()
    return load(str(MEMORY/name/"training.f64"),str(PROBES/(name+".f64")),
                str(MEMORY/name/"rules.f64"),groups)


def inputs():
    parent,plan,_ = tracer_inputs()
    if sha(MEMORY/"summary.json")!=MEMORY_SHA:
        raise ValueError("Saved travel-time memory run changed")
    memory = validate_record(MEMORY/"summary.json")
    if not memory["analysis"]["passed"]:
        raise ValueError("Travel-time memory is unqualified")
    _,f,sp,batch = source_inputs()
    length,area,flow = (f[k] for k in ("length_m","area_m2","flow_m3_s"))
    liquid = sp["porosity"]*length*area
    mass = f["derived_dry_mass_kg"]
    np.testing.assert_allclose(mass,parent["analysis"]["dry_carbon_mass_kg"],rtol=1e-13)
    np.testing.assert_allclose(mass,sp["bulk_density_kg_m3"]*length*area,rtol=1e-13)
    side = parent["analysis"]["native"]["raw"]["side_m"]
    cases = {}
    for size,(_,diameter) in CASES.items():
        width = side*diameter/.0008
        block = round(length/(8*width))
        film = call("film_coefficient",diameter_m=diameter,superficial_velocity_m_s=flow/area,
            viscosity_Pa_s=sp["viscosity_Pa_s"],density_kg_m3=sp["density_kg_m3"],
            molecular_diffusivity_m2_s=sp["molecular_diffusivity_m2_s"],
            porosity=sp["porosity"],interstitial=False,multiplier=1.)
        p = dict(blocks=block,liquid_volume_m3=liquid,dry_mass_kg=mass,flow_m3_s=flow,
            inlet_concentration_kg_m3=f["feed_kg_m3"],particle_diameter_m=diameter,
            particle_density_kg_m3=sp["particle_density_kg_m3"],diffusivity_m2_s=batch["diffusivity_m2_s"],
            film_m_s=film["film_m_s"],capacity_kg_kg=sp["capacity_kg_kg"],
            affinity_m3_kg=sp["affinity_m3_kg"],t_end_s=18000.,startup_s=24.)
        width_used = length/(8*block)
        cases[size] = dict(parameters=p,film=film,reference_section_length_m=width,
            column_section_length_m=width_used,column_section_Pe=(flow/area)/sp["porosity"]*
            width_used/sp["molecular_diffusivity_m2_s"],
            probe_Pe=plan["cases"][f"{size}_primary"]["Pe"])
    return dict(bed_length_m=length,bed_area_m2=area,porosity=sp["porosity"],cases=cases),memory


def solve(name,groups,wall,p):
    model = transport(name,groups)
    _,solver = julia_models()
    fractions = np.array(model.flow)
    areas = model.area if wall else np.tile(fractions[:,None],(1,8))
    r = solver(model.sections,model.flow,areas,**p)
    jl = runtime()
    return {str(f):np.array(v,dtype=float).tolist() if bool(jl.isa(v,jl.AbstractArray)) else float(v)
            for f in jl.fieldnames(jl.typeof(r)) for v in (getattr(r,str(f)),)}


def verify(result,p):
    validate_solution(result,p,dict(dry_inventory_kg=p["dry_mass_kg"],liquid_volume_m3=p["liquid_volume_m3"]))
    area = 6*p["dry_mass_kg"]/(p["particle_density_kg_m3"]*p["particle_diameter_m"])
    np.testing.assert_allclose(result["external_particle_area_m2"],area,rtol=1e-12,atol=0.)
    np.testing.assert_allclose(result["total_film_conductance_m3_s"],area*p["film_m_s"],rtol=1e-12,atol=0.)
    net = result["cumulative_film_uptake_kg"]-result["cumulative_film_release_kg"]
    np.testing.assert_allclose(net,result["adsorbed_inventory_kg"][-1],rtol=1e-8,atol=1e-14)


def inert_trial(name,p,observed):
    model = transport(name,4)
    # An inlet step in the unscaled octet yields the first-passage CDF.
    tau = np.sum(np.array(model.tau))
    test = dict(p,blocks=1,liquid_volume_m3=p["flow_m3_s"]*tau,film_m_s=0.,
                inlet_concentration_kg_m3=1.,t_end_s=6.,dt_s=.0005,startup_s=0.,
                radial_cells=1,seam_history=True)
    r = solve(name,4,True,test)
    verify(r,test)
    np.testing.assert_allclose(r["delay_scale"],1.,rtol=1e-14)
    time = np.array(r["time_s"])
    centres = (time[:-1]+time[1:])/2
    times = np.r_[centres[0]/2,(centres[:-1]+centres[1:])/2]
    cdf = np.array(r["seam_concentration"])[:,1:]
    weights = np.diff(cdf,axis=1,prepend=0.)
    if np.min(weights)<-1e-9 or not np.allclose(cdf[:,-1],1.,rtol=0.,atol=1e-8):
        raise ValueError("Unresolved or nonmonotone inert first passages")
    mean = weights@times
    variance = weights@(times**2)-mean**2
    relative = (np.column_stack((mean,variance))[2:]/np.array(observed)-1)
    maxima = np.max(abs(relative),axis=0)
    if maxima[0]>.05 or maxima[1]>.10:
        raise ValueError(f"The coupled inert model failed at {name}: {maxima}")
    return dict(maximum_mean_error=float(maxima[0]),maximum_variance_error=float(maxima[1]),
        mean_s=mean.tolist(),variance_s2=variance.tolist(),time_step_s=test["dt_s"],
        mass_residual=r["normalized_mass_residual"])


def prepare(base,bed,qualification):
    inert = {}
    for size,setup in bed["cases"].items():
        for step in ("primary","half_step"):
            name = f"{size}_{step}"
            p = setup["parameters"]
            data = qualification["analysis"]["cases"][name]["observed"]
            inert[name] = inert_trial(name,p,data)
            print("Inert transport",name,inert[name]["maximum_variance_error"],flush=True)
            model = transport(name,4)
            setup[step] = dict(wall_hits=np.array(model.hits).astype(int).tolist(),
                contact_fractions=np.array(model.area).tolist(),flow_fractions=np.array(model.flow).tolist(),
                local_mean_delay_s=np.array(model.tau).tolist())
    write_json(base/"input_model.json",bed)
    write_json(base/"inert.json",inert)


def simulate(base,bed):
    predictions = {}
    for size,setup in bed["cases"].items():
        initial = setup["parameters"]
        feed = initial["inlet_concentration_kg_m3"]
        def save(label,step,groups,wall,p):
            r = solve(f"{size}_{step}",groups,wall,p)
            verify(r,p)
            name = f"{size}_{label}.json"
            write_json(base/name,dict(inputs=p,result=r))
            print(name,crossings(r,feed),"runtime s",round(r["wall_s"],2),flush=True)
            return name,r

        runs,refinement = [],[]
        previous = None
        for i,(radial,dt,early) in enumerate(LEVELS):
            p = dict(initial,radial_cells=radial,dt_s=dt,startup_dt_s=early)
            name,result = save(f"wall_{i}","primary",4,True,p)
            runs.append(name)
            if previous is not None:
                latest = full_refinement(previous,result,feed)
                refinement.append(latest)
                print(size,"refinement",latest,flush=True)
                if latest["passed"] and latest["startup"]["passed"]:
                    break
            previous = result
        primary = refinement[-1]
        numerically_qualified = primary["passed"] and primary["startup"]["passed"]
        controls = {key:save(key,step,groups,wall,p)[0] for key,(step,groups,wall) in CONTROLS.items()}
        predictions[size] = dict(selected=runs[-1],refinement_files=runs,refinement=refinement,
                                 numerically_qualified=numerically_qualified,controls=controls)
    return predictions


def score(base,frozen):
    _,source,_,_ = source_inputs()
    if sha(PIXELS)!=source["pixel_csv_sha256"]:
        raise ValueError("Figure 7 marker extraction changed")
    cases = {}
    for size,(series,_) in CASES.items():
        item = frozen["cases"][size]
        cases[size] = {}
        for key,name in {"wall_memory":item["selected"],**item["controls"]}.items():
            record = read_json(base/name)
            r,p = record["result"],record["inputs"]
            m = compare(r,series)
            cases[size][key] = dict(physical_screen_passed=column_passed(m),metrics=m,
                carbon_kg=r["physical_adsorbent_mass_kg"],normalized_mass_residual=r["normalized_mass_residual"],
                normalized_interface_flux_mismatch=r["normalized_interface_flux_mismatch"],
                contact_area_m2=r["external_particle_area_m2"],delay_scale=r["delay_scale"],
                outlet_fraction_60s=float(np.interp(60.,r["time_s"],r["outlet_concentration_kg_m3"])/
                                           p["inlet_concentration_kg_m3"]))
    return cases


def run():
    bed,qualification = inputs()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid4().hex[:8]
    base = OUTPUT/"runs"/stamp
    base.mkdir(parents=True)
    print(base,flush=True)
    files = [HERE/name for name in ("run.py","handoff.jl","model.jl","test_model.jl","protocol.md")]
    files += [JULIA/name for name in ("FilmParticle.jl","ExchangeColumn.jl","SharedContactColumn.jl")]
    files += [HERE.parent/name for name in ("disordered_travel_memory/model.jl","disordered_travel_memory/run.py",
        "disordered_section_exchange/model.jl","disordered_exchange/model.jl","exchange_column/run.py",
        "contact_column/run.py","shared_contact/run.py","alhamed/study.py","alhamed/source.py","film_particle/metrics.py")]
    sources = source_manifest([*files,PARTICLE,SOURCE])
    prepare(base,bed,qualification)
    cases = simulate(base,bed)
    saved = base/"frozen_predictions.json"
    write_json(saved,dict(model_sha256=sha(base/"input_model.json"),inert_sha256=sha(base/"inert.json"),
        curves_sha256={p.name:sha(p) for p in base.glob("d*.json")},cases=cases))
    digest = sha(saved)
    print("Both column curves frozen before Figure 7 scoring",digest,flush=True)
    scored = score(base,read_json(saved))
    summary = base/"summary.json"
    write_json(summary,dict(schema="pore2process.disordered_memory_column.v1",run_id=stamp,
        source_sha256=sources,artifacts_sha256={p.name:sha(p) for p in sorted(base.iterdir())},
        memory_sha256=MEMORY_SHA,particle_batch_sha256=BATCH_SHA,frozen_predictions_sha256=digest,
        figure7=scored,figure7_fitted_parameters=0,new_OpenLB_runs=0,new_DWSIM_runs=0,versions=runtime_versions()))
    write_json(OUTPUT/"latest.json",dict(summary=summary.relative_to(ROOT).as_posix()))
    print(summary,flush=True)


def check(path):
    record = validate_record(path)
    base = path.parent
    saved = base/"frozen_predictions.json"
    frozen = read_json(saved)
    if (sha(saved)!=record["frozen_predictions_sha256"]
            or sha(base/"input_model.json")!=frozen["model_sha256"]
            or sha(base/"inert.json")!=frozen["inert_sha256"]
            or any(sha(base/name)!=digest for name,digest in frozen["curves_sha256"].items())
            or record["figure7"]!=score(base,frozen)):
        raise ValueError("Frozen inputs, curves or Figure 7 scores differ on replay")
    for size,item in frozen["cases"].items():
        files = [*item["refinement_files"],*item["controls"].values()]
        results = [read_json(base/name) for name in files]
        for result in results:
            verify(result["result"],result["inputs"])
        primaries = [r["result"] for r in results[:len(item["refinement_files"])]]
        feed = results[0]["inputs"]["inlet_concentration_kg_m3"]
        if [full_refinement(a,b,feed) for a,b in zip(primaries,primaries[1:])]!=item["refinement"]:
            raise ValueError(f"Refinement record changed for {size}")
    print(dict(sources=len(record["source_sha256"]),artifacts=len(record["artifacts_sha256"])))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else run()
