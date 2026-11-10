"""Keep a small real-simulator reference independently usable by the tests."""
import shutil
from pore2process_demo.io import sha
from pore2process_demo.io import write_json
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read

def main():
    path=ROOT/'results/operating_map/process/latest.json';summary=__import__('pathlib').Path(read(path)['summary']);s=read(summary)
    if s['status']!='passed':raise ValueError('Only accepted results can become a reference.')
    native=read(summary.parent/'native_summary.json');support=read(summary.parent/'support_summary.json')
    c=next(c for c in s['cases'] if c['request']==dict(inlet_Re=1388.,flux_LMH=12.,fixed_multiplier=1.,support='straight'))
    compact_source={k:support[k] for k in ('supports','protocol','fluid_anchor')}
    compact_source['literature_reference']={'fit':{'fixed_impedance_Pa_s_m':support['literature_reference']['fit']['fixed_impedance_Pa_s_m']}}
    value={'status':'passed','summary_path':str(summary),'summary_sha256':sha(summary),
        'native_summary_path':s['native_path'],'native_summary_sha256':s['native_sha256'],
        'native_cases':len(native['cases']),'process_cases':len(s['cases'])+1,'protocol':s['protocol'],
        'hydraulics':s['hydraulics'],'evidence':s['evidence'],'native_for_replay':{k:native[k] for k in (
            'status','geometry','geometry_coefficient_m4','verified_Re_range','accepted_case')},
        'TRT_support_for_replay':compact_source,'actual_case':c,'comparison':s['comparison'],'selections':[v for v in s['selections'] if v['fixed_multiplier']==1 and v['mapping']=='equal_Re' and v['capacity_fraction']==.95 and v['transfer_factor']==1],
        'selection_scope':'Compact subset: fixed factor 1, equal Re, capacity fraction 0.95, transfer factor 1; complete selections remain in the accepted CSV/summary.',
        'verification':read(summary.parent/'verification.json')}
    value['native_for_replay']['cases']={native['accepted_case']:native['cases'][native['accepted_case']]}
    write_json(ROOT/'examples/filtration/operating_map/reference/operating_map_reference.json',value)
    shutil.copy2(summary.parent/'design_map.png',ROOT/'examples/filtration/operating_map/figures/operating_map.png')
    from .report import update
    update(s,summary)
    print('Saved verified operating-map reference, figure and numerical report.')

if __name__=='__main__':main()
