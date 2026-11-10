"""Separate source contracts for Figure 2; legacy Figure 7 extraction is frozen."""
import csv
import json

import numpy as np

from examples.adsorption.alhamed.source import ROOT, PDF, PDF_SHA, evidence as paper_evidence, sha

PIXELS=ROOT/'literature/sources/alhamed2009/figure2_pixels.csv'
PROTOCOL=ROOT/'examples/adsorption/film_particle/protocol.md'


def observations():
    with PIXELS.open() as stream:
        rows=list(csv.DictReader(stream))
    result={}
    for diameter in [1.5,.8,.45,.225]:
        selected=[r for r in rows if float(r['diameter_mm'])==diameter]
        x=np.array([float(r['pixel_x']) for r in selected])
        y=np.array([float(r['pixel_y']) for r in selected])
        result[f'd{diameter*1000:g}']=dict(diameter_m=diameter*.001,role=selected[0]['role'],
            time_s=(x-74)/(589/3000)*60,loading_kg_kg=(30+(425-y)/(416/20))*.001,
            pixel_x=x,pixel_y=y)
    return result


def evidence():
    if sha(PDF)!=PDF_SHA:
        raise ValueError('Paper changed; re-audit source.')
    e=paper_evidence()
    return dict(source_id=e['source_id'],doi=e['doi'],pdf_sha256=PDF_SHA,
        material=e['material'],source_parameters={'capacity_kg_kg':.1115,'affinity_m3_kg':139.,
            'bulk_density_kg_m3':740.,'porosity':.427,'particle_density_kg_m3':740/(1-.427),
            'temperature_K':295.15,'source_temperature_range_K':[294.15,296.15],
            'molecular_diffusivity_m2_s':1e-9,'viscosity_Pa_s':1/1200,'density_kg_m3':1000.},
        figure_2={'location':'Figure 2, journal p765 / PDF p3',
            'raster_name':'page3-image1-Im4.jpg','raster_dimensions':[682,493],
            'raster_sha256':'06bade2a5a80af4f1230980d535d27cef752dc3dc9745750c44e610e1636281c',
            'pixel_origin':'zero-based x rightward, y downward in unchanged embedded raster',
            'axes':{'x0_pixel':74,'x3000min_pixel':663,'y30_mg_g_pixel':425,'y50_mg_g_pixel':9},
            'pixel_csv_sha256':sha(PIXELS),'initial_concentration_kg_m3':.2,
            'fluid_volume_m3':.00025,'dry_mass_kg':.001,'shaker_rpm':200,
            'selected_points':sum(len(x['time_s']) for x in observations().values()),
            'calibration_diameters_m':[.0008,.0015],'holdout_diameters_m':[.00045,.000225],
            'reading_sensitivity_pixels':2.,'reading_sensitivity_time_s':2/(589/3000)*60,
            'reading_sensitivity_loading_kg_kg':2/(416/20)*.001,
            'reading_role':'chosen image-reading sensitivity, not experimental uncertainty',
            'pointwise_measurement_uncertainty':None,
            'sampling':'Section 2.2 reports 1 mL aliquots; complete sampling schedule unavailable.',
            'unobserved':'Uptake below 30 mg/g and unresolved overlapping early markers.'},
        batch_table_3=e['batch_table_3'],figure_7=e['figure_7'],
        equilibrium_factors={'primary_capacity_multiplier':1.,'primary_affinity_multiplier':1.,
            'diagnostic_capacity_factor':'fit only to training Table 3 rows 4, 5, 6; not used silently in forward runs'},
        inventory={'source_bulk_density_kg_m3':740.,'derived_column_dry_mass_kg':e['figure_7']['derived_dry_mass_kg'],
            'mass_role':'density times bed geometry, not a per-run weighed mass',
            'primary_active_fraction':1.,'chosen_active_fraction_scenarios':[.5,.75],
            'scenario_role':'unmeasured accessibility scenarios, not uncertainty bounds or fitted masses'},
        protocol_sha256=sha(PROTOCOL),unavailable=e['unavailable'])


def main():
    p=ROOT/'literature/sources/alhamed2009/film_particle_extraction.json'
    p.write_text(json.dumps(evidence(),indent=2,allow_nan=False)+'\n')
    print(p)


if __name__=='__main__':
    main()
