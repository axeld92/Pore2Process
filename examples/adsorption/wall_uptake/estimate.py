"""Estimate film conductance from captured flux and surviving occupation."""
import numpy as np

from examples.adsorption.contact_column.inputs import CONTACT, MARKERS, rank_bins
from examples.adsorption.shared_contact.run import model as shared_model
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT

from .prepare import requests

LIMITS = dict(step_total_capture=.005,step_total_conductance_relative=.15,
              second_weighted_capture_error=.01,two_period_survival_error=.01)
SHARED = ROOT/"results/adsorption/shared_contact/runs/20261112T121757Z-1a08a56f"
SHARED_SHA = "184b565f8a08823b519a6a23195d74918da7124c6dec60b443682b13fc78fa37"


def uptake(labels, transit, first_hit, n):
    captured = np.isfinite(first_hit)
    if (not np.isfinite(transit).all() or np.any(transit<=0) or np.isinf(first_hit).any()
            or np.any(first_hit[captured]<0) or np.any(first_hit[captured]>transit[captured]+1e-12)):
        raise ValueError("Invalid first-hit time or transit duration")
    counts = np.bincount(labels,minlength=n)
    hit_count = np.bincount(labels,weights=captured,minlength=n)
    total = np.bincount(labels,weights=transit,minlength=n)
    alive = np.where(captured,np.minimum(first_hit,transit),transit)
    occupation = np.bincount(labels,weights=alive,minlength=n)
    if np.any(counts==0) or np.any(occupation<=0):
        raise ValueError("Insufficient tracer occupation in a path class")
    fraction = hit_count/counts
    mean_survival = occupation/total
    return dict(counts=counts,hit_count=hit_count.astype(int),transit_sum_s=total,
                surviving_sum_s=occupation,capture=fraction,occupation_fraction=mean_survival,
                NTU=fraction/mean_survival)


def reduce_case(base,key,spec,labels,mix,w):
    data = np.fromfile(base/(key+".f64"),dtype="<f8").reshape(-1,4)
    old_path = CONTACT/(key+".f64")
    if sha(old_path)!=spec["reference_sha256"]:
        raise ValueError("Archived contact trajectory changed")
    old = np.fromfile(old_path,dtype="<f8").reshape(-1,8)
    np.testing.assert_array_equal(data[:,:2],old[:,:2])
    np.testing.assert_array_equal(np.isfinite(data[:,2:]),old[:,2:4]>0)
    n = len(w)
    first,second = (uptake(labels[:,i],data[:,i],data[:,i+2],n) for i in range(2))
    survival = np.exp(-first["NTU"])
    hit = np.isfinite(data[:,2:])
    measured_survival = float(np.mean(~hit.any(axis=1)))
    predicted_survival = float(w@(survival*(mix@survival)))
    endpoint_only_survival = float(w@((1-second["capture"])*(mix@(1-first["capture"]))))
    errors = dict(second_weighted_capture_error=float(w@abs(1-survival-second["capture"])),
                  two_period_survival_error=abs(predicted_survival-measured_survival))
    return dict(first={k:v.tolist() for k,v in first.items()},second={k:v.tolist() for k,v in second.items()},
        first_total_capture=float(w@first["capture"]),second_reset_total_capture=float(w@second["capture"]),
        first_predicted_total_capture=float(w@(1-survival)),weighted_NTU=float(w@first["NTU"]),
        measured_two_period_survival=measured_survival,predicted_two_period_survival=predicted_survival,
        second_capture_given_first_survival=float(hit[~hit[:,0],1].mean()),
        endpoint_only_survival_at_fixed_exchange=endpoint_only_survival,
        errors=errors,transfer_passed=all(v<=LIMITS[k] for k,v in errors.items()))


def coefficients(ntu,weights,flow_fraction,p,trace):
    area = 6*trace["dry_inventory_kg"]/(p["particle_density_kg_m3"]*p["particle_diameter_m"])
    flow = p["area_m2"]*p["superficial_velocity_m_s"]
    k = flow*flow_fraction*np.array(ntu)*p["periods"]
    contact_area = area*np.array(weights)
    if np.any((contact_area==0)&(k>0)):
        raise ValueError("Inferred film transfer has no assigned interface area")
    per_area = np.divide(k,contact_area,out=np.zeros_like(k),where=contact_area>0)
    base = area*p["film_m_s"]
    return dict(path_film_m_s=per_area.tolist(),path_bed_conductance_m3_s=k.tolist(),
                total_bed_conductance_m3_s=float(k.sum()),correlation_conductance_m3_s=base,
                ratio_to_correlation=float(k.sum()/base),external_area_m2=area,
                zero_transfer_classes=int(np.count_nonzero(k==0)))


def model(base):
    plan = read_json(base/"plan.json")
    if plan["cases"]!=requests(base) or any(sha(ROOT/name)!=d for name,d in plan["sources"].items()):
        raise ValueError("Wall uptake marker request changed")
    if sha(SHARED/"summary.json")!=SHARED_SHA:
        raise ValueError("Shared-carbon reference changed")
    previous = validate_record(SHARED/"summary.json")
    shared = shared_model()
    if shared!=read_json(SHARED/"model.json"):
        raise ValueError("Shared-carbon inputs changed")
    w = np.array(shared["closure"]["flow_fractions"])
    edges = np.r_[0.,np.cumsum(w)]
    marker = validate_record(MARKERS/"summary.json")
    twins = np.fromfile(MARKERS/"advection.f64",dtype="<f8").reshape(-1,2)
    if not np.isfinite(twins).all():
        raise ValueError("Incomplete diffusion-free twins")
    labels = rank_bins(twins,edges)
    seed = marker["requests"]["advection"][-2:]
    sizes = {}
    for size,setup in shared["sizes"].items():
        matrix = np.array(setup["mix"])
        spectra = {}
        for kind in ("primary","half_step"):
            key = size+"_"+kind
            spec = plan["cases"][key]
            if spec["command"][-2:]!=seed:
                raise ValueError("Uptake marker positions differ from their diffusion-free twins")
            spectra[kind] = reduce_case(base,key,spec,labels,matrix,w)
        a,b = (spectra[k] for k in ("primary","half_step"))
        step = dict(step_total_capture=abs(a["first_total_capture"]-b["first_total_capture"]),
                    step_total_conductance_relative=abs(a["weighted_NTU"]/b["weighted_NTU"]-1))
        film = coefficients(a["first"]["NTU"],setup["weights"]["wall"],w,setup["parameters"],setup["trace"])
        sizes[size] = dict(parameters=setup["parameters"],trace=setup["trace"],mix=setup["mix"],
            weights=setup["weights"]["wall"],uptake=spectra,step_changes=step,
            step_passed=all(v<=LIMITS[k] for k,v in step.items()),coefficients=film,
            controls=previous["figure7"][size])
    return dict(closure=shared["closure"],dry_carbon_mass_kg=shared["dry_carbon_mass_kg"],sizes=sizes,
                parent_shared_sha256=SHARED_SHA,limits=LIMITS)
