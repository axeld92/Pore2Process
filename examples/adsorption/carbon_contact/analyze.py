"""Reduce tracked wall encounters to contact shares and fixed-mass scenarios."""
import argparse
from pathlib import Path
import shutil

import numpy as np

from examples.adsorption.exchange_attribution.run import ALLOCATIONS, fractions
from examples.adsorption.exchange_column.closure import closure
from examples.adsorption.exchange_column.run import source_inputs
from examples.adsorption.path_exchange.predict import ranks
from examples.adsorption.porous_face_handoff.run import archived
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

from .geometry import audit, wall_faces
from .prepare import HERE, MARKERS, OUTPUT, specifications

EVENT = np.dtype([("marker_phase","<i4"),("face_id","<u4"),("reflection_m","<f8")])
LIMITS = dict(pair_recovery=.999,uniform_wall_mean_relative=.10,rank_share_step_absolute=.03)


def partition(exposure, inventory):
    """Distribute observed patch inventory once using shares of local time."""
    total = exposure.sum(axis=1)
    observed = total>0
    resolved = (exposure[observed]/total[observed,None]).T@inventory[observed]
    unresolved = float(inventory[~observed].sum())
    np.testing.assert_allclose(resolved.sum()+unresolved,inventory.sum(),rtol=1e-12,atol=0)
    return resolved,unresolved


def patch_weights(by_face, patch, inventory):
    matrix = np.zeros((len(inventory),by_face.shape[1]))
    np.add.at(matrix,patch,by_face)
    return partition(matrix,inventory)


def exposure_result(events, face_rows, event_label, labels, by_marker, wall, mass):
    groups = int(labels.max())+1
    area = wall["smooth_area_m2"]
    by_face = np.zeros((len(area),groups))
    np.add.at(by_face,(face_rows,event_label),events["reflection_m"])
    group_regulator = by_face.sum(axis=0)
    counts = np.bincount(labels.ravel(),minlength=groups)
    grouped_hits = np.bincount(labels.ravel(),weights=by_marker.ravel()>0,minlength=groups)
    count_solid = wall["solid_voxel_count_by_patch"]
    patch, unresolved = patch_weights(by_face,wall["patch"],count_solid/count_solid.sum())
    grain_count = count_solid.reshape(-1,8).sum(axis=1)
    grain, unresolved_grain = patch_weights(by_face,wall["grain"],grain_count/grain_count.sum())
    return dict(group_doses=counts.tolist(), contact_probability=(grouped_hits/counts).tolist(),
                local_time_fraction=(group_regulator/group_regulator.sum()).tolist(),
                local_time_mean_um=(1e6*group_regulator/counts).tolist(),
                patch_carbon_fraction=patch.tolist(),patch_carbon_mass_g=(patch*mass*1000).tolist(),
                patch_unresolved_mass_g=unresolved*mass*1000,
                grain_carbon_fraction=grain.tolist(),grain_unresolved_mass_g=unresolved_grain*mass*1000),by_face


def overlap(by_face, area):
    visited = by_face.sum(axis=1)>0
    fast,slow = by_face[:,0]>0,by_face[:,-1]>0
    other = by_face[:,:-1].sum(axis=1)>0
    af = area/area.sum()
    share, unresolved = partition(by_face,af)
    return dict(area_visited=float(af[visited].sum()), unresolved_area_fraction=unresolved,
                fast_visited_area=float(af[fast].sum()),slow_visited_area=float(af[slow].sum()),
                shared_slow_and_other_area=float(af[slow & other].sum()),
                shared_fast_and_slow_area=float(af[fast & slow].sum()),
                maximum_potentially_slow_exclusive_area=float(af[~other].sum()),
                slow_only_observed_area=float(af[slow & ~other].sum()),
                allocated_face_area_fraction=share.tolist(),
                slow_face_area_allocation_bounds=[float(share[-1]),float(share[-1]+unresolved)])


def preceding_transit(decile, hits, regulator):
    previous = decile[:,0]
    n = np.bincount(previous,minlength=10)
    second_hits = hits[:,1]>0
    count = np.bincount(previous,weights=second_hits,minlength=10)
    exposure = np.bincount(previous,weights=regulator[:,1],minlength=10)
    first_hits = hits[:,0]>0
    return dict(second_contact_probability_by_first_rank=(count/n).tolist(),
                second_exposure_fraction_by_first_rank=(exposure/exposure.sum()).tolist(),
                first_contact_probability=float(first_hits.mean()),
                second_contact_probability=float(second_hits.mean()),
                contact_in_both=float((first_hits&second_hits).mean()),
                second_contact_given_first=float(second_hits[first_hits].mean()),
                second_contact_given_no_first=float(second_hits[~first_hits].mean()))


def analyze_case(base, name, spec, wall, mass, dx, flow, class_edges, no_diffusion=None):
    path = base/(name+".f64")
    data = np.fromfile(path,dtype="<f8").reshape(-1,8)
    reference = ROOT/spec["reference_pairs"]
    if sha(reference)!=spec["reference_sha256"]:
        raise ValueError("Archived transit pairs changed")
    pair_ref = np.fromfile(reference,dtype="<f8").reshape(-1,2)
    np.testing.assert_array_equal(data[:,:2],pair_ref)
    if len(data)!=spec["count"] or not np.isfinite(data[:,2:]).all() or np.any(data[:,2:]<0):
        raise ValueError("Invalid marker history")
    completed = np.isfinite(data[:,:2]).all(axis=1)
    recovery = float(completed.mean())
    events = np.fromfile(path.with_suffix(path.suffix+".faces"),dtype=EVENT)
    marker = abs(events["marker_phase"].astype(np.int64))-1
    phase = (events["marker_phase"]<0).astype(int)
    if len(events)==0 or np.any(marker<0) or np.any(marker>=len(data)) or not np.all(events["reflection_m"]>0):
        raise ValueError("Invalid face events")
    event_code = 2*marker+phase
    totals = np.bincount(event_code,weights=events["reflection_m"],minlength=2*len(data)).reshape(-1,2)
    hits = np.bincount(event_code,minlength=2*len(data)).reshape(-1,2)
    np.testing.assert_allclose(totals,data[:,4:6],rtol=1e-12,atol=1e-18)
    np.testing.assert_array_equal(hits,data[:,2:4])
    order = np.argsort(wall["face_id"])
    ids = wall["face_id"][order]
    face_index = np.searchsorted(ids,events["face_id"])
    if np.any(face_index>=len(ids)) or not np.array_equal(ids[face_index],events["face_id"]):
        raise ValueError("A tracked contact is outside the solid mask")
    observed_rows = order[face_index]
    expected = 2*spec["diffusivity_m2_s"]*len(ids)*dx**2/flow
    reflected = float(totals.sum()/len(data))
    uniform_error = abs(reflected/expected-1)
    valid_event = completed[marker]
    valid_tracks = data[completed]
    rank = np.full((len(data),2),-1)
    rank[completed] = np.column_stack([ranks(p) for p in valid_tracks[:,:2].T])
    decile = np.minimum(9,10*rank//len(valid_tracks))
    class_bins = np.searchsorted(class_edges[1:],rank/len(valid_tracks),side="right")
    contact_rank,face_by_decile = exposure_result(events[valid_event],observed_rows[valid_event],
        decile[marker,phase][valid_event],decile[completed],hits[completed],wall,mass)
    classes,_ = exposure_result(events[valid_event],observed_rows[valid_event],
        class_bins[marker,phase][valid_event],class_bins[completed],hits[completed],wall,mass)
    slow1,slow2 = decile[completed].T==9
    class_overlap = overlap(face_by_decile,wall["smooth_area_m2"])
    patch_encounters = np.zeros((len(wall["solid_voxel_count_by_patch"]),10))
    np.add.at(patch_encounters,wall["patch"],face_by_decile)
    twin = None
    if no_diffusion is not None:
        # Same inlet draw by marker ID. Rank the diffusion-free twin separately
        # so wall hits on the measured trajectory cannot affect its rank label.
        twin_ranks = np.full((len(data),2),-1)
        twin_ranks[completed] = np.column_stack([ranks(p) for p in no_diffusion[completed].T])
        twin_decile = np.minimum(9,10*twin_ranks//len(valid_tracks))
        twin,twin_faces = exposure_result(events[valid_event],observed_rows[valid_event],
            twin_decile[marker,phase][valid_event],twin_decile[completed],hits[completed],wall,mass)
        twin["face_overlap"] = overlap(twin_faces,wall["smooth_area_m2"])
    checks = dict(completed=recovery>=LIMITS["pair_recovery"],
                  uniform_exposure=uniform_error<=LIMITS["uniform_wall_mean_relative"],
                  all_solid_patches_resolved=contact_rank["patch_unresolved_mass_g"]<=1e-12)
    return dict(recovery=recovery,face_events=len(events),uniform_area_oracle=dict(
                    expected_um=1e6*expected,observed_um=1e6*reflected,relative_error=uniform_error),
                exposure=contact_rank,allocation_40_classes=classes,face_overlap=class_overlap,
                grain_patch_overlap=dict(patches=len(patch_encounters),
                    slow_visited=int(np.count_nonzero(patch_encounters[:,-1]>0)),
                    shared_fast_and_slow=int(np.count_nonzero((patch_encounters[:,0]>0)&(patch_encounters[:,-1]>0))),
                    shared_slow_and_other=int(np.count_nonzero((patch_encounters[:,-1]>0)&(patch_encounters[:,:-1].sum(axis=1)>0)))),
                slow_retention=float(slow2[slow1].mean()),advection_ranked_exposure=twin,
                preceding_transit=preceding_transit(decile[completed],hits[completed],totals[completed]),
                external_end_reflections=int(data[:,6].sum()),post_crossing_hits_in_first_geometry=int(data[:,7].sum()),
                checks=checks),face_by_decile


def compare_cases(a,b):
    changes = {k:float(np.max(abs(np.array(a["exposure"][k])-b["exposure"][k])))
               for k in ("local_time_fraction","patch_carbon_fraction","grain_carbon_fraction")}
    if a["advection_ranked_exposure"] is not None and b["advection_ranked_exposure"] is not None:
        x,y = (r["advection_ranked_exposure"]["patch_carbon_fraction"] for r in (a,b))
        changes["advection_ranked_patch_carbon_fraction"] = float(np.max(abs(np.array(x)-y)))
    return dict(max_decile_changes=changes,passed=max(changes.values())<=LIMITS["rank_share_step_absolute"])


def analysis(base):
    plan = read_json(base/"plan.json")
    if plan["cases"]!=specifications(base) or any(sha(ROOT/name)!=digest for name,digest in plan["sources"].items()):
        raise ValueError("Archived marker settings or source changed")
    info, ijk, _ = archived(96)
    wall = wall_faces(info,ijk)
    flow = info["carrier_directional_flow_m3_s"][0]
    _, f, _, _ = source_inputs()
    mass = f["derived_dry_mass_kg"]
    cl = closure()
    shares = np.array(cl["flow_fractions"])
    edges = np.r_[0.,np.cumsum(shares)]
    tracker = read_json(MARKERS/"summary.json")
    file = MARKERS/"advection.f64"
    command = tracker["requests"]["advection"]
    if sha(file)!=tracker["artifacts_sha256"][file.name]:
        raise ValueError("Archived diffusion-free marker twin changed")
    advection = np.fromfile(file,dtype="<f8").reshape(-1,2)
    if not np.isfinite(advection).all():
        raise ValueError("Diffusion-free twin is incomplete")
    reports = {name:analyze_case(base,name,spec,wall,mass,info["dx_m"],flow,edges,
        no_diffusion=advection if spec["command"][-2:]==command[-2:] else None) for name,spec in plan["cases"].items()}
    cases = {name:report for name,(report,_) in reports.items()}
    pooled = {p:overlap(sum(w for k,(_,w) in reports.items() if k.startswith(p+"_")),wall["smooth_area_m2"])
              for p in ("reference","d800","d450")}
    refine = {p:compare_cases(cases[p+"_primary"],cases[p+"_half_step"]) for p in ("reference","d800","d450")}
    replicate = compare_cases(cases["reference_primary"],cases["reference_replicate"])
    indices = np.r_[0.,np.cumsum(shares)]
    fast = np.diff(np.minimum(indices,.1))/shares
    slow = np.diff(np.maximum(indices,.9))/shares
    old = {}
    for label,power in ALLOCATIONS.items():
        a = fractions(cl,power)
        old[label] = dict(fast_fraction=float(fast@a),slow_fraction=float(slow@a))
    return dict(cases=cases,step_refinement=refine,independent_sample=replicate,pooled_face_overlap=pooled,
                previous_carbon_scenarios=old,total_dry_mass_kg=mass,
                all_contact_checks_passed=(all(all(c["checks"].values()) for c in cases.values())
                    and all(v["passed"] for v in refine.values()) and replicate["passed"]),
                wall=dict(staircase_faces=len(wall["face_id"]),staircase_area_m2=len(wall["face_id"])*info["dx_m"]**2,
                          smooth_area_m2=float(sum(wall["smooth_area_m2"])),solid_voxels=int(sum(wall["solid_voxel_count_by_patch"]))),
                geometry_audit=audit(),limits=LIMITS)


def record(base):
    result = analysis(base)
    files = [HERE/name for name in ("prepare.py","analyze.py","test_analyze.py","track_contacts.jl",
                                   "test_track_contacts.jl","geometry.py","protocol.md")]
    files += [HERE.parent/"path_exchange/track.jl",HERE.parent/"path_exchange/predict.py",
              HERE.parent/"exchange_attribution/run.py",HERE.parent/"exchange_column/closure.py"]
    hashes = source_manifest(files)
    for name in hashes:
        target = base/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    path = base/"summary.json"
    write_json(path,dict(schema="pore2process.carbon-contact.v1",run_id=base.name,
        source_sha256=hashes,artifacts_sha256={str(p.relative_to(base)):sha(p) for p in base.rglob("*") if p.is_file()},
        analysis=result,new_OpenLB_runs=0,new_DWSIM_runs=0,new_column_fitted_parameters=0))
    write_json(OUTPUT/"latest.json",dict(summary=str(path.relative_to(ROOT))))
    print("Carbon contact result",path,result["all_contact_checks_passed"],flush=True)


def check(path):
    saved = validate_record(path)
    if analysis(path.parent)!=saved["analysis"]:
        raise ValueError("Carbon contact reduction differs on replay")
    print(dict(sources=len(saved["source_sha256"]),artifacts=len(saved["artifacts_sha256"]),
               contact_checks_passed=saved["analysis"]["all_contact_checks_passed"]))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument("--prepared",type=Path)
    choices.add_argument("--check",type=Path)
    args = parser.parse_args()
    check(args.check.resolve()) if args.check else record(args.prepared.resolve())
