"""Retain published measurements, reading bounds and transfer assumptions."""
from pathlib import Path
from calibration.audit_crossflow_sources import audit
from coupling.run import sha
from hydraulic.run import write_json
from process_filtration.run import ROOT, read

def extract(require_sources=False):
    path=ROOT/'literature/crossflow_extraction.json'
    original=read(path); reviewed=audit(original,require_sources)
    bio_path=ROOT/'literature/chiu2006_bio_extraction.json';bio=read(bio_path)
    readings=reviewed['figure_readings']['chiu2006_fig6']
    points=[]
    for index,(re,point) in enumerate(zip((1388.,2776.,3701.,4626.),readings,strict=True)):
        value=point['critical_flux_LMH']['nominal']
        # Explicit written values supersede reading midpoints only where given.
        if index==0:value=21.4
        if index==2:value=72.5
        low,high=point['critical_flux_LMH']['reading_interval']
        if not low<=value<=high:raise ValueError('Text anchor disagrees with figure reading.')
        points.append({'source_Re':re,'critical_flux_LMH':value,'reading_interval_LMH':[low,high],
            'flux_provenance':'Section 3.2 text anchor, checked against Fig. 6' if index in (0,2) else 'Fig. 6 digitized marker',
            'Re_provenance':'Section 3.2 written Reynolds values; checked against digitized x interval',
            'digitized_Re_interval':point['x']['reading_interval'],'experimental_point_uncertainty':None})
        if not point['x']['reading_interval'][0]<=re<=point['x']['reading_interval'][1]:
            raise ValueError('Written Reynolds number outside figure reading bounds.')
    sources=[]
    for source in original['sources']:
        sources.append({'id':source['id'],'doi':source['doi'],'files':source['files'],
            'role':{'chiu2006':'Primary biological-feed Jcrit(Re) and protocol; no same-condition retention curve.',
                    'chiu2005':'TiO2 comparison only; feed differs. Not pooled with biological measurements.',
                    'garcia2008':'Pressure/energy arithmetic and separate TiO2 benchmark; Q convention unresolved.',
                    'frederic2018':'Cross-sectional resistance model; arbitrary fouling coefficient cannot supply a measured crossflow benefit.'}[source['id']]})
    return {'schema_version':1,'source_extraction_sha256':sha(path),'source_audit':reviewed['source_checks'],
        'biological_feed_extraction':bio,'biological_feed_extraction_sha256':sha(bio_path),
        'sources':sources,'primary_source_doi':'10.1016/j.seppur.2005.11.028','points':points,
        'conditions':original['figure_readings']['chiu2006_fig6']['conditions'],
        'source_rig_orientation':{'value':'vertically mounted','location':'Chiu 2006 Section 2.1'},
        'specimen':'Seven-channel ceramic microfilter, nominal pore size 0.2 micrometre; star-channel family, not reconstructed geometry.',
        'protocol':'TMP steps 30 min before nonlinearity and 15 min after; average of last time-independent and first time-dependent flux. Maximum TMP 1.2 bar.',
        'reported_reproducibility':original['figure_readings']['chiu2006_fig6']['reported_reproducibility'],
        'reading_bounds_scope':'Chosen independent +/-2 pixel marker/tick sensitivity, not experimental uncertainty.',
        'reference_velocity_pair':{'Re':1388.,'velocity_m_s':.6,'location':'Sections 2.1 and 3.2'},
        'nominal_relation':'Piecewise linear interpolation of critical flux versus source Reynolds number; no extrapolation.',
        'shape_relation':'Monotonic brackets without assuming an interpolation shape; exact source nodes use their own reading intervals.',
        'transfer_assumptions':['The synthetic hydraulic case is horizontal with gravity neglected and uniform external permeate pressure; no hydrostatic profile is inferred from the vertically mounted source rig.',
            'The hydraulic fixed-impedance anchor derives from the earlier paired Zhu specimen extraction. Chiu uses a different membrane/feed; this combination is a scenario, not a matched experimental prediction.',
            'Transfer the source critical-flux curve to a synthetic wall/support at unchanged nominal selective-layer scenario.',
            'Compare equal Reynolds number and equal channel velocity as separate mappings; neither identifies the source flow convention.',
            'A constant transfer multiplier and capacity fraction are explicit scenario choices, not measured uncertainty.',
            'Apply source module-average Jcrit as a conservative local-flux screen: maximum J(x) <= scenario capacity; this local use is an assumption.',
            'Water-carrier properties; biological feed rheology unavailable. COD retention remains a separate constant scenario.',
            'Changing support geometry is assumed not to change the source critical-flux curve except through local flux redistribution.'],
        'unavailable':{'sustainable_flux_LMH':None,'fouling_time_s':None,'same_condition_rejection_curve':None,
            'measured_feed_viscosity_Pa_s':None,'pointwise_experimental_uncertainty':None},
        'no_new_article_required_for_this_bounded_milestone':True}

if __name__=='__main__':
    result=extract();write_json(ROOT/'literature/operating_window_evidence.json',result)
    print('Audited four papers; retained four biological-feed critical-flux points and explicit transfer limits.')
