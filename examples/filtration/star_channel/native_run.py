import argparse
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from pore2process_demo.io import sha
from pore2process_demo.io import write_json
from pore2process_demo.bridge import runtime, runtime_versions
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from .build import main as build
from .geometry import polygon, dimensions
from .native import solve, verify_suite

def execute(output):
    p=read(ROOT/'examples/filtration/star_channel/protocol.json'); meta=build(); jl=runtime()
    base=Path(output)/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8])
    base.mkdir(parents=True); write_json(base/'protocol.json',p)
    cases={}
    jobs=[(f'grid{n}',n,{}) for n in p['native']['resolutions']]
    jobs += [('circle',32,{'circle':True}),('half',32,{'gradient':.5}),
             ('double',32,{'gradient':2.}),('reverse',32,{'gradient':-1.}),
             ('domain',32,{'axial':8}),('tau',32,{'tau':.65})]
    for name,n,kw in jobs:
        print('OpenLB: '+name,flush=True); cases[name]=solve(base/name,n,max_steps=p['native']['max_steps'],**kw)
    operating=[]
    for re in (5.,20.,35.):
        operating.append((f'op{re:g}_n64',64,re,{}))
    operating += [('op20_n48',48,20.,{}),('op20_domain',64,20.,{'axial':8}),
                  ('op20_time',64,20.,{'tau':.53}),('op20_cold',48,20.,{'cold':True})]
    from concurrent.futures import ThreadPoolExecutor
    def run_operating(job):
        name,n,re,options=job;options=dict(options);cold=options.pop('cold',False)
        gradient=re/cases[f'grid{n}']['Re_mean']
        print('OpenLB operating range: '+name,flush=True)
        return name,solve(base/name,n,gradient=gradient,tau=options.pop('tau',.55),
            warm=None if cold else base/f'grid{n}/section.csv',warm_scale=gradient,max_steps=p['native']['max_steps'],**options)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for name,value in pool.map(run_operating,operating):cases[name]=value
    references=[]
    vertices=jl.seval('Matrix{Float64}')(polygon())
    for n in (64,128,256):
        v=jl.duct_stokes_reference(vertices,n)
        references.append({'resolution':n,**{str(k):float(getattr(v,str(k))) for k in jl.fieldnames(jl.typeof(v))}})
    checks=verify_suite(cases,references,p)
    current={f.name:sha(f) for f in (ROOT/'examples/filtration/star_channel/openlb').iterdir() if f.is_file()}
    if current!=meta['source_sha256']:raise ValueError('Application source changed during the numerical run; rerun verification.')
    result={'status':'passed','base':str(base.resolve()),'protocol':p,'cases':cases,'accepted_case':'grid64',
            'references':references,'checks':checks,'geometry':dimensions(polygon()),'build':meta,
            'versions':runtime_versions(),
            'source_sha256':{str(f.relative_to(ROOT)):sha(f) for f in (ROOT/'examples/filtration/star_channel/openlb').iterdir() if f.is_file()}}
    result['artifacts_sha256']={str(f.relative_to(base)):sha(f) for f in base.rglob('*') if f.is_file()}
    write_json(base/'summary.json',result); write_json(Path(output)/'latest.json',{'summary':str((base/'summary.json').resolve())})
    print('Native verification passed: '+str(base),flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/star_channel/native')
    execute(parser.parse_args().output)
