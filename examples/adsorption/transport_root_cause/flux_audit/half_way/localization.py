"""Annotate completed candidate peaks using the existing read-only geometry helper."""
import argparse
import json
from pathlib import Path
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from ..findings import peak_statistics
from .run import HERE,OUTPUT,PARENTS


def compute():
    summary=ROOT/json.loads((OUTPUT/'latest.json').read_text())['summary'];s=json.loads(summary.read_text());base=summary.parent
    evidence={str(summary.relative_to(ROOT)):sha(summary)};comparisons={}
    for n in (32,64):
        key=f'porous{n}'
        if key not in s['cases']:continue
        candidate=base/'native'/key;baseline=PARENTS[n]
        records={}
        for label,directory in [('candidate',candidate),('baseline',baseline)]:
            for name in ('peak.csv','peak.json','raw.json','inventory.csv'):
                path=directory/name;evidence[str(path.relative_to(ROOT))]=sha(path)
                if label=='candidate' and s['artifacts_sha256'][str(path.relative_to(base))]!=sha(path):raise ValueError('Candidate peak evidence changed')
            raw=json.loads((directory/'raw.json').read_text())
            records[label]=dict(localization=peak_statistics(directory/'peak.csv',side_m=raw['side_m'],particle_radius_m=.0004),
                                peak_record=json.loads((directory/'peak.json').read_text()))
        old=records['baseline']['peak_record']['maximum_excursion_kg_m3'];new=records['candidate']['peak_record']['maximum_excursion_kg_m3']
        records['candidate_over_baseline_excursion']=new/old
        comparisons[str(n)]=records
    sources=[Path(__file__),HERE.parent/'findings.py',HERE/'run.py']
    return dict(schema='pore2process.halfway-peak-localization.v1',comparisons=comparisons,evidence_sha256=evidence,
                analysis_source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},new_native_runs=0,new_DWSIM_runs=0,
                caveat='Selected peak snapshots describe individual worst-bound events, not simultaneous all-time minima/maxima. Smooth-sphere gaps annotate the prescribed geometry; the actual configured boundary is the voxel midlink. Localization is association rather than proof of causal origin.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute();dest=HERE/'localization.json'
    if a.check:
        if json.loads(dest.read_text())!=result:raise ValueError('Candidate localization changed')
        print('Half-way peak localization replays; no new native/DWSIM runs')
    else:write_json(dest,result);print(dest)


if __name__=='__main__':main()
