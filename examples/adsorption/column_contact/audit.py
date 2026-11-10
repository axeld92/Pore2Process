"""Source, geometry and inventory audit; no experimental uncertainty invented."""
from pore2process_demo.provenance import verified_source_path
import json
import math
import numpy as np
from examples.adsorption.alhamed.source import ROOT, evidence, sha, PDF
from examples.adsorption.alhamed.study import observations

FROZEN_RUN=ROOT/'results/adsorption/film_particle/runs/20261110T083144Z-e33389b9/summary.json'
PROTOCOL=ROOT/'examples/adsorption/column_contact/protocol.md'


def frozen():
    r=json.loads(FROZEN_RUN.read_text())
    for name,digest in r['source_sha256'].items():
        if sha(verified_source_path(name,digest))!=digest:raise ValueError(f'Frozen particle study source changed: {name}')
    return r


def audit():
    old=frozen();e=evidence();f=e['figure_7']
    A=math.pi*.011**2/4;L=.1;Q=43.8e-6/60;eps=.427;bulk=740.;rho=bulk/(1-eps)
    mass=bulk*A*L;c0=.15;q=.1115*139*c0/(1+139*c0)
    capacity=mass*q+eps*A*L*c0
    windows={}
    for label in ('calibration','holdout'):
        t,y=observations(label)
        # No assumption about the unseen tail. Origin is an explicitly declared
        # clean-bed interpolation assumption; even it is not a measured point.
        integral=float(np.trapezoid(1-np.r_[0.,y],np.r_[0.,t]))
        retained=Q*c0*integral
        missing=capacity-retained
        last_deficit=float(1-y[-1])
        # Conditional lower bound: a monotone tail cannot exceed its final
        # observed deficit. A declining tail can require arbitrarily longer.
        extra=max(missing,0)/(Q*c0*last_deficit) if last_deficit>0 else None
        shifts=[]
        for sx in (-1,1):
            for sy in (-1,1):
                tr=np.maximum(0,t+sx*f['reading_sensitivity_time_s'])
                yr=np.clip(y+sy*f['reading_sensitivity_fraction'],0,1)
                shifts.append(float(Q*c0*np.trapezoid(1-np.r_[0.,yr],np.r_[0.,tr])))
        windows[label]=dict(window_end_s=float(t[-1]),last_fraction=float(y[-1]),
            finite_window_retained_kg=retained,
            retained_loading_if_printed_inventory_kg_kg=retained/mass,
            fraction_of_full_equilibrium_inventory=retained/capacity,
            assumed_tail_retention_required_kg=missing,
            minimum_additional_monotone_tail_s=extra,
            reading_sensitivity_retained_kg=[min(shifts),max(shifts)],
            role='trapezoid reading audit, not saturation uptake or a measured bed mass')
    fits=[]
    for row in e['column_table_5']:
        if row['flow_mL_min']==43.8 and row['length_cm']==10 and row['inlet_mg_L']==150:
            fits.append(dict(diameter_mm=row['diameter_mm'],
                fitted_DL_cm2_min=row['fitted_DL_cm2_min'],reported_Pe=row['reported_Pe'],derived_Pe=row['derived_Pe'],
                role='fitted breakthrough-shape descriptors; not independently measured fluid dispersion',
                used_as_model_input=False))
    return dict(source=dict(doi=e['doi'],pdf_sha256=sha(PDF),
            locations=['PDF p2 / journal p764: Tables 1,2 and Section 2.3',
                'PDF p5 / journal p767: Table 5','PDF p6 / journal p768: Figure 7']),
        printed_bulk_density_kg_m3=bulk,printed_bed_voidage=eps,
        derived_particle_envelope_density_kg_m3=rho,derived_bed_mass_kg=mass,
        density_interpretation='Take the printed bulk label literally as dry mass per bed volume; transfer to both size fractions is an assumption.',
        density_method=None,per_column_weighed_mass_kg=None,per_fraction_voidage=None,
        total_pore_volume_m3_kg=.000456,
        implied_skeletal_density_if_all_conventions_compatible_kg_m3=1/(1/rho-.000456),
        skeletal_diagnostic_role='conditional consistency diagnostic; pore-volume and density measurement bases are not documented sufficiently to prove inconsistency',
        alternative_label_diagnostic=dict(if_740_were_particle_density_bulk_kg_m3=740*(1-eps),
            if_740_were_particle_density_bed_mass_kg=740*(1-eps)*A*L,
            forward_use=False,reason='Contradicts printed bulk label and changes the batch particle density; requires requalification, not a silent inventory correction.'),
        geometry=dict(length_m=L,diameter_m=.011,area_m2=A,volume_m3=A*L,
            flow_m3_s=Q,superficial_velocity_m_s=Q/A,interstitial_velocity_m_s=Q/(eps*A),
            empty_bed_contact_time_s=A*L/Q,void_residence_time_s=eps*A*L/Q,
            diameter_to_particle_ratio={'calibration':.011/.0008,'holdout':.011/.00045}),
        inlet_section=dict(source='Section 2.3: 5 cm top glass beads, glass-wool plugs and bottom mesh; downflow',
            bead_section_length_m=.05,bead_voidage=None,
            liquid_residence_upper_bound_s=A*.05/Q,
            bound_role='geometric upper bound with bead section entirely liquid; excludes unreported tubing/plugs',
            full_50cm_tube_liquid_upper_bound_s=A*.5/Q,
            tubing_holdup_m3=None,tracer_RTD=None),
        equilibrium=dict(source_qmax_kg_kg=.1115,source_b_m3_kg=139.,
            loading_at_feed_kg_kg=q,full_solid_inventory_kg=mass*q,
            full_liquid_inventory_kg=eps*A*L*c0,full_equilibrium_inventory_kg=capacity,
            capacity_area_time_s=capacity/(Q*c0),
            role='asymptotic mass-balance target, not an experimental saturation observation'),
        observed_windows=windows,paper_dispersion=fits,
        prior_batch_study=str(FROZEN_RUN.relative_to(ROOT)),
        prior_batch_source_sha256=old['source_sha256'],
        batch_diffusivity_m2_s=old['batch_fit']['diffusivity_m2_s'],
        batch_film_ref_m_s=old['batch_fit']['batch_film_ref_m_s'],
        applied_previous_column_film_multiplier=False,protocol_sha256=sha(PROTOCOL),
        missing=['weighed mass and voidage for each packed fraction','inert tracer RTD',
            'unobserved breakthrough tails','experimental pointwise uncertainty','column pressure loss'])
