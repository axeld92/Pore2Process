"""Decode visually audited Alhamed marker coordinates with explicit raster axes."""
import csv

import numpy as np

from examples.adsorption.alhamed.source import PDF, PDF_SHA, evidence
from examples.adsorption.alhamed.study import crossing
from examples.adsorption.figure7_retention.run import integrate
from pore2process_demo.io import read_json, sha
from pore2process_demo.paths import ROOT

PIXELS = ROOT/"literature/sources/alhamed2009/figure45_pixels.csv"
EXTRACTION = PIXELS.with_name("figure45_extraction.json")
SERIES = {5:[(10,23.3),(15,23.3),(20,23.3),(25,23.3)],4:[(25,23.3),(25,43.8),(25,62.5)]}


def case_name(depth,flow):
    return f"q{round(10*flow)}_l{depth}"


def operating_cases():
    if sha(PDF)!=PDF_SHA:
        raise ValueError("Alhamed source PDF changed")
    table = evidence()["column_table_5"]
    cases = {}
    for series in SERIES.values():
        for L,Q in series:
            row, = (r for r in table if (r["length_cm"],r["flow_mL_min"],r["diameter_mm"],r["inlet_mg_L"])
                    == (L,Q,.8,150))
            cases[case_name(L,Q)] = dict(bed_cm=L,flow_mL_min=Q,table5_row=row["row"])
    return cases


def markers():
    metadata = read_json(EXTRACTION)
    if sha(PIXELS)!=metadata["pixel_csv_sha256"] or sha(PDF)!=PDF_SHA:
        raise ValueError("Figure 4–5 digitization source changed")
    with PIXELS.open() as file:
        rows = list(csv.DictReader(file))
    data = {}
    for figure,series in SERIES.items():
        axis = metadata["figures"][str(figure)]
        ts = 400*60/(axis["x400min"]-axis["x0"])
        cs = 1/(axis["y0"]-axis["y1"])
        for L,Q in series:
            subset = [r for r in rows if (int(r["figure"]),int(r["bed_cm"]),float(r["flow_mL_min"]))==(figure,L,Q)]
            x,y = (np.array([int(r[k]) for r in subset]) for k in ("pixel_x","pixel_y"))
            time = (x-axis["x0"])*ts
            if len(time)<10 or time[0]<=0 or not np.all(np.diff(time)>0):
                raise ValueError("Unordered or insufficient published markers")
            key = f"f{figure}_{case_name(L,Q)}"
            data[key] = dict(figure=figure,case=case_name(L,Q),time=time,
                             fraction=np.clip((axis["y0"]-y)*cs,0.,1.),
                             dt=metadata["reading_sensitivity_pixels"]*ts,
                             dc=metadata["reading_sensitivity_pixels"]*cs)
    return data


def observed_crossing(t,c,threshold):
    return crossing(t,c,threshold) if c[0]<=threshold else None


def observation(data,feed_rate,carbon,common_time):
    t,c = data["time"],data["fraction"]
    age,conc = np.r_[0.,t],np.r_[0.,c]
    # Compare actual marker windows. The clean-bed origin supplies an integration
    # convention, and does not resolve a crossing before the first marker.
    loading = integrate(age,conc,feed_rate,t[-1])/carbon
    readings = []
    crossings = {".05":[],".50":[]}
    for sx in (-1,1):
        for sy in (-1,1):
            times = np.maximum(0.,t+sx*data["dt"])
            frac = np.clip(c+sy*data["dc"],0.,1.)
            readings.append(float(np.trapezoid(1-np.r_[0.,frac],np.r_[0.,times])*feed_rate/carbon))
            for key in crossings:
                crossings[key].append(observed_crossing(times,frac,float(key)))
    bounds = {k:None if None in v else [min(v),max(v)] for k,v in crossings.items()}
    return dict(marker_count=len(t),end_s=float(t[-1]),last_outlet_fraction=float(c[-1]),
                t05_s=observed_crossing(t,c,.05),t50_s=observed_crossing(t,c,.50),
                retained_kg_kg=loading,reading_retention_kg_kg=[min(readings),max(readings)],
                reading_t05_s=bounds[".05"],reading_t50_s=bounds[".50"],
                common_dose_retained_kg_kg=integrate(age,conc,feed_rate,common_time)/carbon)
