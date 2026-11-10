"""Published facts, conflicts and scenarios stay separate; no pooled specimens."""
import json
import math
from pathlib import Path
from coupling.run import sha
from hydraulic.build import ROOT
from hydraulic.run import write_json

PDF=ROOT/'literature/inbox/Meneghini2017Phenol/MENEGHINI.pdf'
DIGEST='84a6d3c1f196ab0a89e2e83012051fa31164ee744de20f577056227c92296993'

def evidence():
    if PDF.exists() and sha(PDF)!=DIGEST:raise ValueError('Author-source PDF changed; re-audit extraction.')
    rows=[]
    def datum(name,value,unit,si,si_unit,place,role='reported'):
        rows.append(dict(name=name,value=value,unit=unit,SI_value=si,SI_unit=si_unit,
            source='Meneghini2017',location=place,role=role,pointwise_uncertainty=None))
    datum('temperature',25,'degC',298.15,'K','Table 2, p23; Table 4, p26')
    datum('Langmuir_capacity',41.69,'mg/g',.04169,'kg/kg dry solid','Table 2, p23','reported fitted batch parameter')
    datum('Langmuir_affinity',.094,'L/mg',94.,'m3/kg solute','Table 2, p23','reported fitted batch parameter')
    datum('effective_intraparticle_diffusivity_at_50mg_L',6.5973e-8,'cm2/s',6.5973e-12,'m2/s','Table 4, p25','reported fitted parameter; not independent prediction')
    datum('molecular_diffusivity',9.1e-6,'cm2/s',9.1e-10,'m2/s','Table 4, p25','reported model input')
    datum('particle_diameter',.23,'cm',.0023,'m','Table 4 continuation, p26')
    datum('bed_length',10.,'cm',.1,'m','Table 4 continuation, p26')
    datum('column_diameter',2.04,'cm',.0204,'m','Table 4 continuation, p26')
    datum('bed_porosity',.44,'1',.44,'1','Table 4 continuation, p26')
    datum('dry_adsorbent_mass',20.,'g',.020,'kg','Table 4, p25')
    datum('density_as_printed',574.33,'g/L',574.33,'kg/m3','Table 4, p25','density interpretation unresolved')
    datum('column_flow',20.,'mL/min',20e-6/60,'m3/s','Table 4, p25')
    V=.1*math.pi*.0204**2/4
    return {'status':'literature-informed starting parameters; experimental column validation not passed or claimed',
        'sources':[{'id':'Meneghini2017','type':'author undergraduate dissertation','title':'Simulação numérica e ensaios experimentais da remoção de fenol em coluna de leito fixo utilizando carvão ativado como adsorvente',
            'material':'Carbomafra coconut-shell activated carbon','source_url':'https://rd.uffs.edu.br/handle/prefix/1681',
            'pdf_url':'https://rd.uffs.edu.br/bitstream/prefix/1681/1/MENEGHINI.pdf','sha256':DIGEST,
            'local_path':str(PDF.relative_to(ROOT)),'available_locally':PDF.exists(),'visual_audit_pages':[23,25,26]},
            {'id':'Meneghini2021','doi':'10.1590/S1413-415220180046','relationship':'Related journal article; not used as independent validation of the author-source parameters.',
             'access':'publisher HTML reviewed; PDF not retrieved','source_url':'https://www.scielo.br/j/esa/a/k97n87Dwkz5RLKMKCW3Rd3M/?lang=pt'},
            {'id':'Alhamed2009','doi':'10.1016/j.jhazmat.2009.05.002','material':'date-stone activated carbon',
             'access':'publisher abstract only','role':'Candidate particle-size validation series; no numbers borrowed into coconut-shell model.'},
            {'id':'Iheanacho2021','doi':'10.1016/j.sajce.2021.02.003','material':'corn-cob activated carbon',
             'access':'publisher abstract only','role':'Alternative particle-size series; separate specimen, not pooled.'}],
        'measurements_and_fitted_parameters':rows,
        'source_consistency':{'bed_volume_m3':V,'bulk_density_from_mass_and_dimensions_kg_m3':.020/V,
            'particle_density_if_mass_dimensions_porosity_are_used_kg_m3':.020/(V*(1-.44)),
            'mass_if_printed_density_is_particle_density_kg':574.33*(1-.44)*V,
            'printed_density_not_used_as_an_unambiguous_particle_density':True},
        'conflicts':['Table 2 affinity is L/mg; Table 4 prints L/g for a similar number. Use only the unambiguous Table 2 units.',
            'Journal Table 2 has 44.15 mg/g and 0.07 L/mg, differing from author Table 2 and journal Table 5; do not merge fits.',
            'Printed density, mass, porosity and dimensions cannot all satisfy one particle-envelope inventory definition.',
            'Figure 7 early effluent and error bars require digitization and an initial/dead-volume audit before a clean-bed holdout comparison.',
            'Source radial-diffusion model is more detailed than the initial intraparticle LDF approximation.'],
        'uncertainty_status':'Figure 7 has visible error bars; numerical intervals and their definition have not yet been digitized/audited. The related journal methodology reports duplicate tests and 95% confidence; no pointwise standard deviation is inferred.',
        'unavailable':{'experimental_pressure_drop':None,'pointwise_measurement_uncertainty':None,
            'matched_particle_size_effect_for_this_carbon':None,'validated_bed_dispersion_m2_s':None},
        'scenario_assumptions':{'synthetic_geometry':'periodic FCC sphere bed, porosity 0.6; not reconstruction of carbon granules',
            'particle_envelope_dry_density_kg_m3':1000.,'bed_length_m':1.,'bed_area_m2':.01,
            'breakthrough_Cout_Cin':.05,'particle_diameters_m':[.000575,.00115,.0023],
            'axial_dispersion_m2_s':0.,'film_resistance':'not separately resolved in this initial LDF model',
            'dispersion_sensitivity_factors':[.1,1.],'dispersion_sensitivity_definition':'Dax = chosen factor * interstitial velocity * particle diameter; scenario only, not fitted uncertainty.',
            'kinetic_transfer':'60 Deff/dp^2 assumes spherical intraparticle control and unchanged material Deff.',
            'DWSIM_role':'Actual water-carrier pump; adsorption, solute inventory and transient capacity remain external Julia calculations.'}}

if __name__=='__main__':
    e=evidence();write_json(ROOT/'literature/adsorption_evidence.json',e)
    print(json.dumps({'sources':len(e['sources']),'retained_entries':len(e['measurements_and_fitted_parameters']),'consistency':e['source_consistency']},indent=2))
