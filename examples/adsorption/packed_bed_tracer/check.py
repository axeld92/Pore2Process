"""Verify saved native evidence and replay Julia, without launching OpenLB."""
import json
from pathlib import Path
import shutil
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.runtime.openlb import verify_vendor
from .analyze import analyze
from .run import relative,translate


def hashes(root,manifest):
    root=Path(root).resolve()
    for name,digest in manifest.items():
        p=(root/name).resolve()
        if not p.is_relative_to(root) or sha(p)!=digest:raise ValueError(f'Changed or unsafe artifact: {name}')


def check(summary):
    path=ROOT/summary;base=path.parent;s=json.loads(path.read_text())
    hashes(ROOT,s['source_sha256']);hashes(base/'source_snapshot',s['source_sha256']);hashes(base,s['artifacts_sha256'])
    cases={}
    for key,expected in s['cases'].items():
        p=base/'native'/key;request=json.loads((p/'request.json').read_text());spec=request['spec']
        assert json.loads((p/'execution.json').read_text())['returncode']==0
        source=spec['flow_source']
        if source:
            assert sha(ROOT/source['raw_path'])==source['raw_sha256']
            assert sha(spec['field'])==source['field_sha256']
            assert sha(ROOT/'results/adsorption/bed/summary.json')==source['source_summary_sha256']
        cases[key]=analyze(p)
        assert cases[key]==expected
        final=np.genfromtxt(p/'final.csv',delimiter=',',names=True)
        r=expected['raw'];assert len(final)==r['nx']*r['n']**2
        assert len(set(zip(final['i'],final['j'],final['k'])))==len(final)
        mask=(final['i']>0)&(final['i']<r['nx']-1)
        fluid=mask&(final['material']!=2);walls=mask&(final['material']==2)
        vol=r['dx_m']**3
        np.testing.assert_allclose(np.sum(final['concentration_kg_m3'][fluid])*vol+
            np.sum(final['concentration_kg_m3'][walls])*vol,expected['final_inventory_kg'],rtol=1e-10,atol=1e-20)
        np.testing.assert_allclose(np.sum(fluid)*vol,r['fluid_volume_m3'],rtol=1e-10)
    refinement=dict(grid=relative(cases['bed32'],cases['bed24']),time=relative(cases['time32'],cases['bed32']))
    assert refinement==s['refinements']
    assert all(all(c['checks'].values()) for c in cases.values())==s['native_checks_passed']
    assert all(v['mean_s']<=.05 and v['variance_s2']<=.1 for v in refinement.values())==s['numerical_refinement_passed']
    # Replay to a verification-only temporary directory; never overwrite native
    # outputs or their already-hashed Julia predictions.
    import tempfile
    with tempfile.TemporaryDirectory(dir=ROOT/'.local') as tmp:
        tmp=Path(tmp);(tmp/'native').symlink_to(base/'native',target_is_directory=True)
        contact,predictions=translate(tmp,cases)
    assert contact==s['contact']
    for key,actual in predictions.items():
        expected=json.loads((base/f'julia_{key}.json').read_text())
        assert actual==expected
        assert actual['metrics']==s['macro_comparison'][key]
    assert s['new_native_scalar_runs']==8 and s['new_native_flow_runs']==1 and s['new_DWSIM_runs']==0
    vendor=verify_vendor();assert vendor==s['vendor']
    ref={k:s[k] for k in ('schema','run_id','cases','refinements','native_checks_passed','numerical_refinement_passed',
        'contact','macro_comparison','versions','new_native_scalar_runs','new_native_flow_runs','new_DWSIM_runs')}
    ref['summary_path']=str(path.relative_to(ROOT))
    dest=ROOT/'examples/adsorption/packed_bed_tracer/reference';dest.mkdir(exist_ok=True)
    write_json(dest/'tracer_reference.json',ref)
    report=dict(saved_native_cases_checked=len(cases),new_native_solves_in_checker=0,
        Julia_predictions_replayed=len(predictions),artifacts_checked=len(s['artifacts_sha256']),
        source_files_checked=len(s['source_sha256']),native_checks_passed=s['native_checks_passed'],
        numerical_refinement_passed=s['numerical_refinement_passed'],macro_comparison=s['macro_comparison'])
    write_json(base/'verification.json',report);print(json.dumps(report,indent=2))
    return report

def main():
    p=json.loads((ROOT/'results/adsorption/packed_bed_tracer/latest.json').read_text());check(p['summary'])

if __name__=='__main__':main()
