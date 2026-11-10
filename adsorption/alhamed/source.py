"""Audited source values; never merge these with the coconut-shell dataset."""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PDF = ROOT / 'literature/inbox/Alhamed2009Adsorption/Alhamed2009Adsorption.pdf'
PDF_SHA = 'c0d9eb61560f643b5b890b1926f1b0f2bc279386dc6bd882cb640ae3f0061b22'
PIXELS = ROOT / 'literature/alhamed2009_fig7_pixels.csv'
RASTER_SHA = '09f8685d2898e06efbb530f2d19e1a5632a75931911d9988f5a460cbe10ef2c8'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def points():
    """Zero-based pixel centres in the unchanged 674 x 532 embedded Im9 raster."""
    with PIXELS.open() as stream:
        rows = list(csv.DictReader(stream))
    out = {}
    for series in ('calibration', 'holdout'):
        selected = [r for r in rows if r['series'] == series]
        x = np.array([float(r['pixel_x']) for r in selected])
        y = np.array([float(r['pixel_y']) for r in selected])
        # x: 0 min at 83 px; 400 min at 606 px. y: 0 at 471; 1 at 9.
        out[series] = {'diameter_m': float(selected[0]['diameter_mm']) * 1e-3,
                       'time_s': (x - 83) / (523 / 400) * 60,
                       'fraction': (471 - y) / 462,
                       'pixel_x': x, 'pixel_y': y}
    return out


def evidence():
    if sha(PDF) != PDF_SHA:
        raise ValueError('Alhamed source changed; re-audit the PDF.')
    batch = [
        # dp mm, C mg/L, dose g, second-order qe mg/g, k2 g/(mg min), R2, MSE, h
        (.225, 200, 1., 48.49, .02588, .99995, .0163, 60.85),
        (.45, 200, 1., 48.76, .00448, .99999, .0045, 10.65),
        (.45, 400, 1., 90.3, .00207, .9998, .2004, 16.88),
        (.8, 200, 1., 48.95, .00209, .99987, .0406, 5.01),
        (1.5, 200, 1., 48.75, .00074, .99997, .0085, 1.76),
        (1.5, 400, 1., 90.4, .00022, .99870, 1.111, 1.80),
        (.45, 200, .5, 85.4, .00080, .99984, .135, 5.83),
    ]
    raw_columns = [
        # Q mL/min,L cm,C mg/L,dp mm,FLUB,DL cm2/min,Vmin cm,R2,Pe
        (23.9,5,150,1.47,1.,2096,147.8,.983,.140),
        (42.6,5,150,1.47,1.,19867,29.1,.919,.026),
        (43.8,10,150,.8,.923,1601,1445.1,.997,.674),
        (62.5,15,150,.8,.985,2515,2479.1,.997,.918),
        (43.8,10,100,.8,.977,1260,2203.8,.998,.856),
        (141.5,25,150,.8,.954,8711,3910.7,.996,1.),
        (43.8,10,50,.8,.969,1214,3518.6,.997,.888),
        (43.8,15,150,.8,.947,1615,2687.2,.997,1.002),
        (88.1,25,150,.8,.896,3913,4616.2,.996,1.386),
        (23.3,10,150,.8,.923,428,2114.2,.999,1.339),
        (62.5,20,150,.8,.896,2272,4103.8,.998,1.355),
        (62.5,25,150,.8,.923,2541,4792.2,.998,1.515),
        (43.8,10,150,.45,.799,171,2369.2,1.,6.309),
        (43.8,20,150,.8,.861,1150,4233.,.999,1.876),
        (43.8,25,150,.8,.817,1025,5979.9,.997,2.631),
        (23.3,15,150,.8,.783,285,3721.2,.998,3.021),
        (23.3,20,150,.8,.783,405,4829.2,.999,2.832),
        (23.3,25,150,.8,.637,218,6758.9,.999,6.594),
    ]
    columns = []
    area = math.pi * .011**2 / 4
    for i, (Q,L,C,dp,flub,D,V,R2,Pe) in enumerate(raw_columns, 1):
        u = Q*1e-6/60/area
        columns.append(dict(row=i,flow_mL_min=Q,length_cm=L,inlet_mg_L=C,
            diameter_mm=dp,FLUB=flub,fitted_DL_cm2_min=D,fitted_Vmin_cm=V,
            fitted_R2=R2,reported_Pe=Pe,derived_Pe=(u/.427)*(L*.01)/(D*1e-4/60),
            location='Table 5, journal p767 / PDF p5',
            role='reported operating conditions and fitted column parameters',
            pointwise_uncertainty=None))
    measurements = []
    for name,value,unit,si,si_unit,location,role in [
        ('bulk_density',.74,'g/cm3',740.,'kg dry solid/m3 packed bed','Table 1, p764','reported characteristic'),
        ('bed_porosity',.427,'1',.427,'1','Table 1, p764','reported characteristic'),
        ('Langmuir_capacity',111.5,'mg/g',.1115,'kg/kg dry solid','Table 1, p764','reported fit from reference 5'),
        ('Langmuir_affinity',.139,'L/mg',139.,'m3/kg solute','Table 1, p764','reported fit from reference 5'),
        ('Freundlich_k',24.6,'(mg/g)(L/mg)^(1/n)',None,None,'Table 1, p764','reported fit; not used in this Langmuir test'),
        ('Freundlich_n',2.74,'1',2.74,'1','Table 1, p764','reported fit from reference 5'),
        ('column_diameter',11.,'mm',.011,'m','Section 2.3 / Table 2, p764','reported apparatus'),
        ('temperature',22.,'degC',295.15,'K','Sections 2.2 and 2.3, p764','reported condition; +/- 1 degC'),
        ('liquid_viscosity',.5,'g/(cm min)',.5*.001/.01/60,'Pa s','Table 2, p764','reported model input'),
        ('molecular_diffusivity',1e-5,'cm2/s',1e-9,'m2/s','Table 2, p764','reported model input'),
        ('water_density',1.,'g/cm3',1000.,'kg/m3','Table 2, p764','reported model input'),
    ]:
        measurements.append(dict(name=name,value=value,unit=unit,SI_value=si,
            SI_unit=si_unit,location=location,role=role,uncertainty=None,
            temperature_range_K=[294.15,296.15] if name=='temperature' else None))
    return dict(source_id='Alhamed2009',doi='10.1016/j.jhazmat.2009.05.002',
        title='Adsorption kinetics and performance of packed bed adsorber for phenol removal using activated carbon from dates’ stones',
        material='ZnCl2-activated date-stone carbon; separate from coconut-shell carbon',
        pdf_path=str(PDF.relative_to(ROOT)),pdf_sha256=PDF_SHA,visual_audit_pdf_pages=[2,3,4,5,6,7],
        measurements_and_fits=measurements,
        batch_table_3=[dict(experiment=i,diameter_mm=d,inlet_mg_L=C,dose_g=m,
            fitted_qe_mg_g=q,fitted_k2_g_mg_min=k,fitted_R2=r,fitted_MSE=e,
            fitted_initial_rate_mg_g_min=h,solution_volume_L=.25,
            location='Table 3, p765',role='published batch fits, not raw observations',
            pointwise_uncertainty=None) for i,(d,C,m,q,k,r,e,h) in enumerate(batch,1)],
        column_table_5=columns,
        figure_7=dict(location='Figure 7, p768 / PDF p6',
            calibration_diameter_m=.0008,holdout_diameter_m=.00045,
            length_m=.1,flow_m3_s=43.8e-6/60,feed_kg_m3=.150,
            area_m2=area,derived_dry_mass_kg=740.*area*.1,
            digitization='manual centres of open experimental markers, not solid model lines',
            pixel_origin='zero-based x rightward, y downward in original embedded Im9 JPEG',
            raster_dimensions=[674,532],raster_sha256=RASTER_SHA,
            axes={'x0_pixel':83,'x400min_pixel':606,'y0_pixel':471,'y1_pixel':9},
            pixel_csv_sha256=sha(PIXELS),reading_sensitivity_pixels=2.,
            reading_sensitivity_time_s=2/(523/400)*60,
            reading_sensitivity_fraction=2/462,
            reading_sensitivity_role='chosen reading scenario; not experimental uncertainty',
            pointwise_uncertainty=None),
        cautions=[
            'Figure 6 caption says L=20 cm, but matching 50/100 mg/L Table 5 rows say L=10 cm; exclude that series from this test.',
            'Largest particle is called 1.47 mm in methods/Table 5 and 1.5 mm in batch Table 3/abstract; retain each as printed.',
            'Table 5 includes R2=0.983 and 0.919 despite the blanket high-R2 claim; do not discard these rows.',
            'Table 2 omits 15 cm, whereas Table 5 contains 15 cm cases; retain actual row conditions.',
            'Equation 10 visibly prints Pe/4, not sqrt(Pe/4); volume-area convention also needs auditing. No silent correction.',
            'Fitted DL changes with size, bed height and concentration; it cannot be assumed independently measured hydrodynamic dispersion.',
            'Table 4 segment diffusivities depend strongly on size and segment and have nonzero intercepts; not one transferable Deff.',
            'No measured pressure drop or pointwise Figure 7 uncertainty is supplied.',
            'Particle-scale pore volume/BET do not define a connected resolved specimen geometry.',
        ],unavailable={'pressure_drop_Pa':None,'pointwise_measurement_uncertainty':None,
                       'experimental_3D_connected_geometry':None,'pH':None})


def main():
    path=ROOT/'literature/alhamed2009_extraction.json'
    path.write_text(json.dumps(evidence(),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print(path)


if __name__=='__main__':
    main()
