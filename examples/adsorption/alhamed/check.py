"""Replay the paper comparison, retaining a failed physical screen as a failure."""
from pore2process_demo.provenance import verified_source_path
import json
from pathlib import Path

import numpy as np

from examples.adsorption.bridge import column
from .source import ROOT, evidence, sha
from .study import compare, write, source_equation_audit


def check(path):
    path=Path(path);base=path.parent;s=json.loads(path.read_text())
    if s['evidence']!=evidence():
        raise ValueError('Current source extraction differs from the recorded comparison.')
    for name,digest in s['source_sha256'].items():
        if sha(verified_source_path(name,digest))!=digest or sha(base/'source_snapshot'/name)!=digest:
            raise ValueError(f'Changed comparison source: {name}')
    for name,digest in s['artifacts_sha256'].items():
        if sha(base/name)!=digest:
            raise ValueError(f'Changed comparison artifact: {name}')
    replays=0
    for file in base.glob('*.json'):
        c=json.loads(file.read_text())
        if 'result' not in c or 'inputs' not in c:
            continue
        actual=column(**c['inputs']);old=c['result']
        for key,value in old.items():
            if key=='total_s':
                continue
            if value is None:
                assert actual[key] is None
            else:
                np.testing.assert_allclose(actual[key],value,rtol=1e-9,atol=1e-13,err_msg=f'{file.name}: {key}')
        if actual['normalized_mass_residual']>1e-8:
            raise ValueError('Replayed conservation failed.')
        serie='calibration' if file.name.startswith('calibration') else 'holdout'
        m=compare(actual,serie)
        assert m==c['metrics']
        replays+=1
    assert replays==10
    assert source_equation_audit()==s['source_equation_10_audit']
    calibration=next(iter(s['metrics']['calibration'].values()))
    holdout=next(iter(s['metrics']['holdout'].values()))
    calibrated=calibration['RMSE_screen_passed'] and calibration['t50_screen_passed']
    numerical=all(v<.02 for series,c in s['refinement'].items() for r in c.values()
                  for k,v in r.items() if k=='t50' or k=='t05' and series=='holdout')
    numerical &= all(r['normalized_mass_residual']<=1e-8
                     for c in s['refinement'].values() for r in c.values())
    early=s['early_calibration_t05_refinement']
    numerical &= early['space_relative']<.02 and early['time_relative']<.02
    numerical &= all(early[k]['normalized_mass_residual']<=1e-8 for k in ['base','space','time'])
    early_results={k:json.loads((base/f'calibration-early-{k}.json').read_text()) for k in ['base','space','time']}
    for k,c in early_results.items():
        assert early[k]['t05_s']==c['metrics']['t05']['predicted_s']
        assert early[k]['normalized_mass_residual']==c['result']['normalized_mass_residual']
    assert early['space_relative']==abs(early['space']['t05_s']/early['base']['t05_s']-1)
    assert early['time_relative']==abs(early['time']['t05_s']/early['base']['t05_s']-1)
    passed=numerical and calibrated and holdout['RMSE_screen_passed'] and holdout['t50_screen_passed'] and holdout['t05_screen_passed'] is not False
    assert s['calibration_passed']==calibrated
    assert s['numerical_refinement_passed']==numerical
    assert s['physical_comparison_passed']==passed
    return dict(status='passed_evidence_and_replay_checks',Julia_columns_replayed=replays,
                source_files_checked=len(s['source_sha256']),artifacts_checked=len(s['artifacts_sha256']),
                numerical_refinement_passed=numerical,physical_comparison_passed=passed,
                experimental_hydraulic_validation=False)


def main():
    path=Path(json.loads((ROOT/'results/adsorption/alhamed2009/latest.json').read_text())['summary'])
    result=check(path);write(path.parent/'verification.json',result)
    s=json.loads(path.read_text())
    reference={k:s[k] for k in ['status','source_pdf_sha256','fit','metrics','refinement',
        'numerical_refinement_passed','calibration_passed','physical_comparison_passed',
        'source_equation_10_audit','versions','limitations','early_calibration_t05_refinement','refinement_scope']}
    reference['summary_path']=str(path);reference['verification']=result
    reference['cases']=[json.loads((path.parent/name).read_text()) for name in ['calibration-1.json','holdout-1.json','holdout-2.json']]
    for case in reference['cases']:
        case['result']={k:v for k,v in case['result'].items() if not isinstance(v,list) and k!='total_s'}
    write(ROOT/'examples/adsorption/alhamed/reference/alhamed2009_reference.json',reference)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
