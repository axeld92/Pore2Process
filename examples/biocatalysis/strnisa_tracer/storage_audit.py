"""Nonreactive accessible-liquid and diffusion-storage bounds.

The related Lipozyme 435 specimen is not Strnisa's Novozym 435 lot.
These are deliberately separate moment scenarios, not a fitted tracer law.
"""
import numpy as np
from pore2process_demo.io import read_json,sha
from pore2process_demo.paths import ROOT

EVIDENCE=ROOT/"examples/biocatalysis/pentyl_propanoate_transport/particle/evidence.json"


def parameters():
    e=read_json(EVIDENCE)["lorenzo_2025"]
    return dict(source=e["article"],specimen=e["specimen"],
                evidence_path=str(EVIDENCE.relative_to(ROOT)),evidence_sha256=sha(EVIDENCE),
                porosity_scenarios={"dry_sorption":e["particle_porosity"]["dry_nitrogen_sorption"]["value"],
                                    "aqueous_tracers":e["particle_porosity"]["aqueous_tracers"]["value"]},
                borrowed_tortuosity=e["effectiveness_factor"]["adopted_tortuosity"],
                note="Independent related support specimen, not glucose access measurements in Strnisa. Partition 1, no film, uniform connected spherical liquid pores. Moment-only sensitivity.")


def check(plan,reactor,ranking,params):
    """Mean storage follows volume balance regardless of its access rate.

    For linear diffusion in a sphere with prescribed surface concentration,
    the dimensionless mean-concentration transfer function is
    3*(z*coth(z)-1)/z**2 = 1 - z**2/15 + ... . In a plug-flow,
    independent storage approximation its delay variance is
    2*(accessible liquid volume / Q)*(R**2/(15*D_pore)).
    """
    p=plan["protocol"]
    geometry=plan["reactors"][reactor]
    vtotal=float(np.prod(geometry["extent_m"]))
    vsolid=vtotal*(1-geometry["geometry"]["continuous_porosity"])
    q=p["fluid"]["flow_ul_min"]*1e-9/60
    radius=p["geometry"]["diameter_um"]*1e-6/2
    dpore=p["fluid"]["glucose_diffusivity_m2_s"]/params["borrowed_tortuosity"]
    tau=radius**2/(15*dpore)
    observed=ranking["incremental_observed"]
    bed=ranking["incremental_predicted"]
    dm=observed["mean_s"]-bed["mean_s"]
    dv=observed["variance_s2"]-bed["variance_s2"]
    cases={}
    for name,porosity in params["porosity_scenarios"].items():
        added=porosity*vsolid/q
        mean=bed["mean_s"]+added
        variance=bed["variance_s2"]+2*added*tau
        cases[name]=dict(porosity=porosity,additional_mean_s=added,additional_variance_s2=2*added*tau,
                         total_mean_s=mean,total_variance_s2=variance,
                         incremental_mean_relative=mean/observed["mean_s"]-1,
                         incremental_variance_relative=variance/observed["variance_s2"]-1)
    return dict(channel_volume_ul=vtotal*1e9,bead_envelope_volume_ul=vsolid*1e9,
                apparent_experimental_liquid_volume_ul=observed["mean_s"]*q*1e9,
                unexplained_mean_s=dm,unexplained_variance_s2=dv,unexplained_liquid_volume_ul=dm*q*1e9,
                implied_accessible_fraction_of_bead_envelope=dm*q/vsolid,
                diffusion_time_s=tau,assumed_pore_concentration_diffusivity_m2_s=dpore,
                cases=cases,can_positive_storage_explain_mean=dm>=0)
