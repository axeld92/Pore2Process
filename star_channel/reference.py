"""Retain a compact, verified numerical reference without full generated runs."""
import argparse
from pathlib import Path
import shutil

from coupling.run import sha
from hydraulic.run import write_json
from process_filtration.run import ROOT,read
from .check import verify

def retain(path):
    path=Path(path).resolve();s=read(path);verification=verify(path,upstream=False)
    native=read(path.parent/'native_summary.json');source=read(path.parent/'source_summary.json')
    actual=next(c for c in s['cases'] if c['request']==dict(support='straight',flux_LMH=12.,
        rejection=.65,viscosity_multiplier=1.,fixed_multiplier=1.,inlet_Re=20.))
    n={'accepted_case':native['accepted_case'],'geometry':native['geometry'],
       'cases':{k:native['cases'][k] for k in ('grid48','grid64')}}
    support={k:source[k] for k in ('protocol','fluid_anchor','supports','literature_reference')}
    support['supports']=[{k:v[k] for k in ('name','bulk_radius_m','body_resistance_Pa_s_m3')} for v in source['supports']]
    support['literature_reference']={'fit':source['literature_reference']['fit']}
    support['native_cases']={k:{'raw':v['raw']} for k,v in source['native_cases'].items()}
    result={'status':'passed','summary_path':str(path),'summary_sha256':sha(path),
        'timestamp_utc':s['timestamp_utc'],'protocol':s['protocol'],'hydraulics':s['hydraulics'],
        'versions':s['versions'],'verification':verification,'native_verification':s['native_verification'],
        'support_verification':s['support_verification'],'native_for_replay':n,
        'TRT_support_for_replay':support,'actual_case':actual,'comparison':s['comparison'],
        'native_summary_sha256':s['native_sha256'],'support_summary_sha256':s['source_sha256'],
        'scope':'Compact actual numerical reference; full native fields and saved DWSIM cases stay in the uniquely named generated run.'}
    write_json(ROOT/'docs/star_channel_reference.json',result)
    shutil.copy2(path.parent/'study.png',ROOT/'docs/star_channel_study.png')
    shutil.copy2(path.parent/'process'/actual['id']/'flowsheet.png',ROOT/'docs/star_channel_flowsheet.png')
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('summary',nargs='?',type=Path)
    args=parser.parse_args();path=args.summary or Path(read(ROOT/'results/star_channel/process/latest.json')['summary'])
    retain(path);print('Retained verified reference and scientific figures in docs/.')
