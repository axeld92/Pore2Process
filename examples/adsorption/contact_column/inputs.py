"""Build both rank-based allocations from the archived solid-wall events."""
import numpy as np

from examples.adsorption.carbon_contact.analyze import EVENT, partition
from examples.adsorption.carbon_contact.geometry import wall_faces
from examples.adsorption.carbon_contact.prepare import ATTRIBUTION, ATTRIBUTION_SHA, MARKERS
from examples.adsorption.exchange_attribution.run import receiving
from examples.adsorption.exchange_column.closure import closure
from examples.adsorption.exchange_column.run import CASES, source_inputs
from examples.adsorption.path_exchange.predict import ranks
from examples.adsorption.porous_face_handoff.run import archived
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT

CONTACT = ROOT/"results/adsorption/carbon_contact/runs/20261112T101823Z-a87ddcbc"
CONTACT_SHA = "e7d1d50d38612b649ee213e86bda0371fc8cfc0f5d212dfecb813d170955a760"


def rank_bins(times, edges):
    n = len(times)
    limits = np.rint(n*edges).astype(int)
    labels = np.searchsorted(limits[1:],np.column_stack([ranks(t) for t in times.T]),side="right")
    for column in labels.T:
        np.testing.assert_array_equal(np.bincount(column,minlength=len(limits)-1),np.diff(limits))
    return labels


def shares_from_events(events, wall, marker, phase, labels, count):
    order = np.argsort(wall["face_id"])
    ids = wall["face_id"][order]
    index = np.searchsorted(ids,events["face_id"])
    if np.any(index>=len(ids)) or not np.array_equal(ids[index],events["face_id"]):
        raise ValueError("Wall contact does not belong to the saved solid geometry")
    count_by_patch = wall["solid_voxel_count_by_patch"]
    exposure = np.zeros((len(count_by_patch),count))
    np.add.at(exposure,(wall["patch"][order[index]],labels[marker,phase]),events["reflection_m"])
    shares, unresolved = partition(exposure,count_by_patch/count_by_patch.sum())
    if unresolved>0:
        raise ValueError("Contact data leave unresolved octant mass")
    np.testing.assert_allclose(shares.sum(),1.,rtol=1e-12,atol=0)
    return shares


def case_allocations(name, wall, edges, twins, contact, expected_counts):
    trajectory_file = CONTACT/(name+".f64")
    times = np.fromfile(trajectory_file,dtype="<f8").reshape(-1,8)[:,:2]
    if len(times)!=len(twins) or not np.isfinite(times).all():
        raise ValueError("Contact and diffusion-free marker pairs must be complete and aligned")
    events = np.fromfile(trajectory_file.with_suffix(".f64.faces"),dtype=EVENT)
    marker = abs(events["marker_phase"].astype(np.int64))-1
    phase = (events["marker_phase"]<0).astype(int)
    if np.any(marker<0) or np.any(marker>=len(times)) or not np.all(events["reflection_m"]>0):
        raise ValueError("Invalid archived wall events")
    advective = rank_bins(twins,edges)
    realized = rank_bins(times,edges)
    if expected_counts is not None:
        observed = np.zeros_like(expected_counts)
        np.add.at(observed,(realized[:,0],realized[:,1]),1)
        np.testing.assert_array_equal(observed,expected_counts)
    n = len(edges)-1
    advective_mass = shares_from_events(events,wall,marker,phase,advective,n)
    realized_mass = shares_from_events(events,wall,marker,phase,realized,n)
    np.testing.assert_allclose(realized_mass,contact["allocation_40_classes"]["patch_carbon_fraction"],rtol=1e-12,atol=1e-15)
    return dict(advection_rank=advective_mass.tolist(),realized_rank=realized_mass.tolist())


def model():
    for base,digest in ((CONTACT,CONTACT_SHA),(ATTRIBUTION,ATTRIBUTION_SHA)):
        if sha(base/"summary.json")!=digest:
            raise ValueError("Qualified parent cohort changed")
    contact = validate_record(CONTACT/"summary.json")["analysis"]
    prior = validate_record(ATTRIBUTION/"summary.json")
    if not contact["all_contact_checks_passed"]:
        raise ValueError("Carbon contact checks failed")
    marker_archive = validate_record(MARKERS/"summary.json")
    twins = np.fromfile(MARKERS/"advection.f64",dtype="<f8").reshape(-1,2)
    if not np.isfinite(twins).all():
        raise ValueError("The paired advection reference is incomplete")
    contact_plan = read_json(CONTACT/"plan.json")["cases"]
    seed_and_count = marker_archive["requests"]["advection"][-2:]
    _, source, _, _ = source_inputs()
    mass = source["derived_dry_mass_kg"]
    np.testing.assert_allclose(contact["total_dry_mass_kg"],mass,rtol=1e-12,atol=0)
    geometry, ijk, _ = archived(96)
    wall = wall_faces(geometry,ijk)
    cl = closure()
    attribution = read_json(ATTRIBUTION/"model.json")
    if cl!=attribution["closure"]:
        raise ValueError("Exchange class definition changed")
    w = np.array(cl["flow_fractions"])
    edges = np.r_[0.,np.cumsum(w)]
    last = np.diff(np.maximum(edges,.9))/w
    first = np.diff(np.minimum(edges,.1))/w
    sizes = {}
    for size in CASES:
        counts = np.array(attribution["marker_laws"][size]["primary"]["counts"])
        matrix = receiving(counts)
        np.testing.assert_array_equal(matrix,attribution["mix"]["matched_Pe"][size])
        primary, half_step = size+"_primary", size+"_half_step"
        if any(contact_plan[k]["command"][-2:]!=seed_and_count for k in (primary,half_step)):
            raise ValueError("Inlet seeds or marker counts do not match")
        full = case_allocations(primary,wall,edges,twins,contact["cases"][primary],counts)
        half = case_allocations(half_step,wall,edges,twins,contact["cases"][half_step],None)
        summaries = {}
        for rank, fractions in full.items():
            frac = np.array(fractions)
            diff = frac-np.array(half[rank])
            summaries[rank] = dict(fastest_tenth_carbon_fraction=float(first@frac),
                slowest_tenth_carbon_fraction=float(last@frac),zero_contact_classes=int(np.count_nonzero(frac==0)),
                half_step_total_variation=float(.5*np.sum(abs(diff))),
                half_step_slowest_tenth_change=float(last@diff))
        old = prior["analysis"][size]["scenarios"]["void"]["matched_Pe"]
        sizes[size] = dict(mix=matrix.tolist(),allocations=full,allocation_summary=summaries,
                           half_step_allocations=half,void_proportional_control=old)
    return dict(carbon_contact_sha256=CONTACT_SHA,exchange_attribution_sha256=ATTRIBUTION_SHA,
                dry_carbon_mass_kg=mass,closure=cl,sizes=sizes)
