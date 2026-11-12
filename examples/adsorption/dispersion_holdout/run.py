"""Fit the one-cell variance and validate its frozen two-cell prediction."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from examples.adsorption.porous_outlet_tracer.n384 import RESULTS as N384_RESULTS, inputs as outlet_inputs
from examples.adsorption.porous_outlet_tracer.run import PULSE_S, endpoint, history
from examples.adsorption.porous_tracer_refinement.mesh import coarse_faces
from examples.adsorption.porous_tracer_refinement.time_half import CENTER
from examples.adsorption.stationary_curved_wall import reference as ref
from examples.adsorption.transport_accuracy.run import validate_record
from pore2process_demo.io import read_json, sha, write_json
from pore2process_demo.paths import ROOT
from pore2process_demo.provenance import source_manifest

HERE = Path(__file__).parent
OUTPUT = ROOT/'results/adsorption/dispersion_holdout'
PARENT = N384_RESULTS/'20261111T211006Z-4ab30fe9/summary.json'
PARENT_SHA = 'c45bc2e0bfecd2df1b08f4d3b574e4e886eb3cad14d2b567576e5f10b34428d5'
DT = 8e-5
RMSE_LIMIT = .03


def inputs():
    if sha(PARENT)!=PARENT_SHA:
        raise ValueError('The qualified one-cell pulse changed')
    parent = validate_record(PARENT)
    if not parent['analysis']['short_segment_qualified']:
        raise ValueError('The archived one-cell outlet check failed')
    old,runtime,info,coarse = outlet_inputs()
    return parent,runtime,info,coarse,old['steps']['288']


def repeat_mesh(path,info,coarse,*,write):
    source = PARENT.parent/'compact'
    original = {name:np.fromfile(source/name,dtype=dtype).reshape(-1,3) for name,dtype in
                [('before.i32','<i4'),('after.i32','<i4'),('incoming.f64','<f8'),('outgoing.f64','<f8')]}
    n = len(coarse)
    face = coarse_faces(CENTER.parent/'N96')
    l,r,axis = (face[k]-1 for k in ('left.i32','right.i32','axis.u8'))
    seam = (axis==0)&(coarse[l,0]==95)
    out,enter,q = l[seam],r[seam],face['balanced.f64'][seam]
    np.testing.assert_array_equal(coarse[out,1:],coarse[enter,1:])
    np.testing.assert_array_equal(coarse[enter,0],0)
    if np.any(q<=0):
        raise ValueError('Unexpected reverse flow at outer faces')
    np.testing.assert_array_equal(original['outgoing.f64'][out,0],q)
    np.testing.assert_array_equal(original['incoming.f64'][enter,0],q)
    np.testing.assert_array_equal(original['after.i32'][out,0],0)
    np.testing.assert_array_equal(original['before.i32'][enter,0],0)
    arrays = {name:np.tile(a,(2,1)) for name,a in original.items()}
    for name in ('before.i32','after.i32'):
        a = arrays[name]
        a[n:] += n*(a[n:]>0)
    arrays['after.i32'][out,0] = n+enter+1
    arrays['before.i32'][n+enter,0] = out+1
    if write:
        path.mkdir()
    for name,a in arrays.items():
        if write:
            a.tofile(path/name)
        else:
            np.testing.assert_array_equal(np.fromfile(path/name,dtype=a.dtype).reshape(-1,3),a)
    np.testing.assert_allclose(q.sum(),info['carrier_directional_flow_m3_s'][0],rtol=1e-12,atol=0)
    return dict(periods=2,coarse_cells=2*n,fine_cells=2*n*27,dx_m=info['dx_m']/3,
                length_m=192*info['dx_m'],volume_m3=2*n*info['dx_m']**3,
                flow_m3_s=float(q.sum()),coarse_connected_seam_faces=len(q))


def execute(cmd,log,timeout):
    start = perf_counter()
    with log.open('w') as out:
        result = subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,timeout=timeout)
    if result.returncode:
        raise RuntimeError(f'Command failed. See {log}')
    return dict(returncode=result.returncode,wall_s=perf_counter()-start)


def requests(base,parent,runtime,mesh,short):
    variance = parent['analysis']['N384']['variance_s2']
    model = [*runtime,str(HERE/'reference.jl'),*map(str,[base,mesh['length_m']/2,mesh['volume_m3']/2,
              mesh['flow_m3_s'],variance,PULSE_S,DT,short,2*short])]
    dose = mesh['flow_m3_s']*PULSE_S
    holdout = [runtime[0],'--threads=10',*runtime[1:],str(HERE/'solve.jl'),
               *map(str,[base/'compact',base/'holdout',base/'model_2.csv',mesh['dx_m'],ref.D,DT,dose,2*short])]
    return dict(model=model,holdout=holdout,steps=2*short,dose_kg=dose)


def normalized(h,dose):
    return h['recovered_kg']/dose


def analyze_model(base,parent,mesh,request):
    dose = request['dose_kg']
    v = mesh['flow_m3_s']*mesh['length_m']/mesh['volume_m3']
    d = float((base/'coefficient.txt').read_text())
    data = [history(base/f'model_{n}.csv') for n in (1,2)]
    coarse = [history(base/f'model_{n}_coarse.csv') for n in (1,2)]
    tau = mesh['volume_m3']/mesh['flow_m3_s']
    metrics = {}
    for n,h,lo in zip((1,2),data,coarse,strict=True):
        np.testing.assert_array_equal(h['time_s'],lo['time_s'])
        max_cdf = float(max(abs(normalized(h,dose)-normalized(lo,dose))))
        analytical_mean = n*tau/2+PULSE_S/2
        pe = v*mesh['length_m']*n/(2*d)
        analytical_variance = 2*d*d/v**4*(pe+np.expm1(-pe))+PULSE_S**2/12
        stats = endpoint(h,n*tau/2)
        if (max_cdf>1e-4 or abs(stats['mean_exit_s']/analytical_mean-1)>1e-3
                or abs(stats['variance_s2']/analytical_variance-1)>1e-3):
            raise ValueError('Reference curve does not resolve its analytical moments')
        metrics[str(n)] = dict(stats=stats,max_refinement_cdf_change=max_cdf,
            analytical_mean_s=analytical_mean,analytical_variance_s2=float(analytical_variance),Peclet=float(pe))
    expected = parent['analysis']['N384']['variance_s2']
    np.testing.assert_allclose(metrics['1']['analytical_variance_s2'],expected,rtol=1e-12,atol=0)
    archived = history(PARENT.parent/'N384/history.csv')
    np.testing.assert_array_equal(data[0]['time_s'],archived['time_s'])
    short_rmse = float(np.sqrt(np.mean((normalized(data[0],dose)-normalized(archived,dose))**2)))
    return dict(D_m2_s=d,velocity_m_s=v,short_variance_target_s2=expected,reference=metrics,
                short_cdf_rmse=short_rmse)


def analyze_holdout(base,mesh,request):
    path = base/'holdout'
    h = history(path/'history.csv')
    pred = history(base/'model_2.csv')
    n = len(h)
    dose = request['dose_kg']
    np.testing.assert_array_equal(h['time_s'],pred['time_s'][:n])
    np.testing.assert_allclose(h['injected_kg'],mesh['flow_m3_s']*np.minimum(h['time_s'],PULSE_S),
                               rtol=1e-12,atol=0)
    fraction = normalized(h,dose)
    model_fraction = normalized(pred,dose)[:n]
    lower = np.sqrt(np.cumsum((fraction-model_fraction)**2)/(request['steps']+1))
    status = (path/'status.txt').read_text().strip()
    mass = float(max(abs(h['liquid_kg']+h['recovered_kg']-h['injected_kg']))/dose)
    lo,hi = float(min(h['min_kg_m3'])),float(max(h['max_kg_m3']))
    if not np.isfinite(h.view(np.float64)).all() or mass>1e-8 or lo< -1e-8 or hi>1+1e-8:
        raise ValueError('Holdout failed mass, finiteness, or concentration checks')
    first = np.flatnonzero((lower>RMSE_LIMIT)&(fraction>=.1))
    if status=='profile_rejected':
        if not len(first) or first[0]!=n-1:
            raise ValueError('Wrong predeclared profile rejection event')
    elif status!='complete' or n!=len(pred):
        raise ValueError('Incomplete holdout without verified rejection')
    early = {str(f):float(np.interp(f,fraction,h['time_s'])) if fraction[-1]>=f else None for f in (.05,.1)}
    prediction = {str(f):float(np.interp(f,normalized(pred,dose),pred['time_s'])) for f in (.05,.1)}
    result = dict(status=status,physical_time_s=float(h['time_s'][-1]),samples=n,planned_samples=len(pred),
                  recovered_fraction=float(fraction[-1]),model_fraction_at_last_sample=float(model_fraction[-1]),
                  profile_rmse_lower_bound=float(lower[-1]),profile_limit=RMSE_LIMIT,
                  early_arrival_s=early,predicted_early_arrival_s=prediction,max_mass_residual=mass,
                  concentration_range_kg_m3=[lo,hi],qualified=False,full_window_moments=None)
    if status=='complete':
        measured = endpoint(h,mesh['volume_m3']/mesh['flow_m3_s'])
        modeled = endpoint(pred,mesh['volume_m3']/mesh['flow_m3_s'])
        mean = abs(modeled['mean_exit_s']/measured['mean_exit_s']-1)
        variance = abs(modeled['variance_s2']/measured['variance_s2']-1)
        result['full_window_moments'] = dict(measured=measured,modeled=modeled,
                                            mean_relative_error=mean,variance_relative_error=variance)
        result['qualified'] = bool(lower[-1]<=RMSE_LIMIT and mean<=.05 and variance<=.10 and measured['tail_passed'])
    return result


def run():
    parent,runtime,info,coarse,short = inputs()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base = OUTPUT/'runs'/stamp
    base.mkdir(parents=True)
    files = [HERE/k for k in ('run.py','reference.jl','solve.jl','protocol.md')]
    files += [ROOT/'julia/adsorption'/k for k in ('ColumnBoundaryReference.jl','BrickOpenTracer.jl','OpenBoundedFaceTransport.jl')]
    hashes = source_manifest(files)
    for name in hashes:
        target = base/'source_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target)
    print('Dispersion holdout',base,flush=True)
    mesh = repeat_mesh(base/'compact',info,coarse,write=True)
    request = requests(base,parent,runtime,mesh,short)
    record = dict(schema='pore2process.dispersion-holdout.v1',run_id=stamp,source_sha256=hashes,
                  parent_summary_sha256=PARENT_SHA,mesh=mesh,request=request,new_openlb_runs=0,new_bed_runs=0)
    write_json(base/'plan.json',record)
    macro_run = execute(request['model'],base/'reference.log',600)
    model = analyze_model(base,parent,mesh,request)
    frozen = dict(model=model,files={name:sha(base/name) for name in ('coefficient.txt','model_1.csv','model_2.csv')})
    write_json(base/'frozen_prediction.json',frozen)
    frozen_hash = sha(base/'frozen_prediction.json')
    print('Frozen one-cell fit and two-cell prediction',model,flush=True)
    (base/'holdout').mkdir()
    holdout_run = execute(request['holdout'],base/'holdout/stdout.log',28800)
    write_json(base/'execution.json',dict(reference=macro_run,holdout=holdout_run))
    if sha(base/'frozen_prediction.json')!=frozen_hash or any(sha(base/k)!=v for k,v in frozen['files'].items()):
        raise ValueError('Frozen prediction changed after holdout launch')
    record.update(frozen_prediction_sha256=frozen_hash,analysis=analyze_holdout(base,mesh,request))
    record['artifacts_sha256'] = {str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    path = base/'summary.json'
    write_json(path,record)
    write_json(OUTPUT/'latest.json',dict(summary=str(path.relative_to(ROOT))))
    print('Completed two-cell holdout',path,record['analysis'],flush=True)


def check(path):
    path = path.resolve()
    record = validate_record(path)
    parent,runtime,info,coarse,short = inputs()
    base = path.parent
    mesh = repeat_mesh(base/'compact',info,coarse,write=False)
    request = requests(base,parent,runtime,mesh,short)
    frozen = read_json(base/'frozen_prediction.json')
    if (record['parent_summary_sha256']!=PARENT_SHA or mesh!=record['mesh'] or request!=record['request']
            or sha(base/'frozen_prediction.json')!=record['frozen_prediction_sha256']
            or frozen['model']!=analyze_model(base,parent,mesh,request)
            or any(sha(base/k)!=v for k,v in frozen['files'].items())
            or record['analysis']!=analyze_holdout(base,mesh,request)):
        raise ValueError('Holdout differs on replay')
    result = dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),sources=len(record['source_sha256']),
                  artifacts=len(record['artifacts_sha256']),qualified=record['analysis']['qualified'])
    write_json(OUTPUT/'verification.json',result)
    print(result)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path)
    args = parser.parse_args()
    check(args.check) if args.check else run()
