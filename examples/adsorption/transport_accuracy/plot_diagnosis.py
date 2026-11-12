"""Plot the bounded diagnosis; no solver, input or acceptance change."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pore2process_demo.io import sha, write_json
from .run import HERE, OUTPUT, validate_record


def main():
    ref = HERE/'reference/diagnosis_reference.json'
    s = json.loads(ref.read_text())
    from pore2process_demo.paths import ROOT
    native = ROOT/s['source_summary']
    if sha(native) != s['source_summary_sha256']:
        raise ValueError('Changed diagnosis summary.')
    validate_record(native)
    keys = ['reference','no_advection','half_dt','tenth_dt','grid48_matched_dt']
    records = {'reference':s['original_case'], **s['cases']}
    labels = ['N32\noriginal dt','N32\nno flow','N32\nhalf dt','N32\ntenth dt','N48\nhalf dt']
    values = [max(.25-records[k]['minimum_concentration_kg_m3'],
                  records[k]['maximum_concentration_kg_m3']-.75,0.) for k in keys]
    fig,ax = plt.subplots(figsize=(8.6,4.5),layout='constrained')
    ax.bar(np.arange(5),values,color=['#777777','#1b7837','#b2182b','#b2182b','#d98a22'],width=.65)
    ax.axhline(.02,color='.2',ls='--',label='Chosen engineering allowance: 0.02 kg/m³')
    for i,value in enumerate(values):
        ax.text(i,value+.004,f'{value:.3f}',ha='center',fontsize=10)
    ax.set_xticks(np.arange(5),labels)
    ax.set_ylim(0,max(values)*1.23)
    ax.set_ylabel('Largest excursion outside initial range (kg/m³)')
    ax.set_title('Smooth tracer in the periodic porous bed: bounded diagnosis')
    ax.grid(axis='y',alpha=.2)
    ax.set_axisbelow(True)
    ax.legend(fontsize=9,loc='upper left')
    out = HERE/'figures/bounded_diagnosis.png'
    fig.savefig(out,dpi=180)
    plt.close(fig)
    write_json(OUTPUT/'diagnosis_figure_sources.json',dict(reference_sha256=sha(ref),
        native_summary_sha256=sha(native),plotting_source_sha256=sha(Path(__file__)),
        image_sha256=sha(out),visual_inspection='pending'))
    print('Saved',out)


if __name__=='__main__':
    main()
