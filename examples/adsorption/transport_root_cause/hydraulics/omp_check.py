"""Cross-check stock OpenMP CPU_SISD against the same serial native case."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
from time import perf_counter
from uuid import uuid4

import numpy as np

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.provenance import source_manifest
from pore2process_demo.runtime.openlb import OPENLB,verify_vendor
from examples.adsorption.bed import analyze
from .run import HERE,OUTPUT,latest_parent,native_checks


def build_omp(base):
    build=ROOT/'.local/adsorption/transport_root_cause_flow_omp'
    build.mkdir(parents=True,exist_ok=True)
    core=build/'core';core.mkdir(exist_ok=True)
    for source in (HERE/'openlb').iterdir():
        if source.is_file():shutil.copy2(source,build/source.name)
    config=build/'config.mk'
    config.write_text(config.read_text().replace('PARALLEL_MODE := OFF','PARALLEL_MODE := OMP')+'\nOMPFLAGS := -fopenmp\n')
    compiler='/opt/homebrew/bin/g++-15'
    flags=['-O2','-Wall','-std=c++20','-pthread','-fPIC','-DPLATFORM_CPU_SISD',
        '-DPARALLEL_MODE_OMP','-fopenmp','-DOLB_VERSION="1.9r0"',
        f'-I{OPENLB/"src"}','-I/opt/homebrew/include']
    sources=('communication/mpiManager.cpp','communication/ompManager.cpp',
        'core/olbInit.cpp','core/expr.cpp','io/ostreamManager.cpp')
    commands=[];objects=[];start=perf_counter()
    with (base/'build.log').open('w') as log:
        for relative in sources:
            obj=core/(relative.replace('/','_')+'.o');objects.append(obj)
            command=[compiler,*flags,'-c',str(OPENLB/'src'/relative),'-o',str(obj)]
            commands.append(command)
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        command=['ar','rcs',str(core/'libolbcore.a'),*map(str,objects)]
        commands.append(command);subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
        (build/'capillary3d.o').unlink(missing_ok=True)
        command=['make','-j2','GIT_VERSION=1.9r0',f'OLB_ROOT={OPENLB}',f'PROJECT_ROOT={ROOT}','onlysample']
        commands.append(command)
        subprocess.run(command,cwd=build,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    metadata=dict(commands=commands,core_flags=flags,build_s=perf_counter()-start,
        config=(build/'config.mk').read_text(),compiler=subprocess.check_output([compiler,'--version'],text=True),
        vendor=verify_vendor())
    write_json(base/'build.json',metadata)
    binary=base/'binary/capillary3d';binary.parent.mkdir()
    shutil.copy2(build/'capillary3d',binary)
    return binary,metadata


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--serial-binary',type=Path,required=True)
    parser.add_argument('--threads',type=int,default=4)
    args=parser.parse_args()
    if not 1<=args.threads<=8:raise ValueError('Use one to eight native OpenMP threads.')
    hydraulic=latest_parent()
    parent=json.loads((hydraulic/'summary.json').read_text())
    carrier=parent['inputs']
    accepted=parent['cases']['d800']['native']['32']['accepted']
    native=hydraulic/accepted['directory']
    raw=accepted['analysis']['raw']
    rid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
    base=OUTPUT/'omp_checks'/rid;base.mkdir(parents=True)
    files=[Path(__file__),HERE/'run.py',HERE/'protocol.md',*list((HERE/'openlb').iterdir()),
        ROOT/'examples/adsorption/transport_root_cause/flux_audit/population_audit.h',
        ROOT/'examples/adsorption/bed.py']
    sources=source_manifest(files)
    for name in sources:
        dest=base/'source_snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    print('OpenMP cross-check',base,flush=True)
    binary,build=build_omp(base)
    cases={}
    for label,b,threads in [('serial',args.serial_binary.resolve(),1),('omp',binary,args.threads)]:
        dest=base/label;dest.mkdir()
        command=[str(b),'32',repr(raw['gradient_Pa_m']),repr(raw['tau_even']),'.0008',str((native/'field.csv').resolve()),'1.']
        environment=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE')
        write_json(dest/'request.json',dict(command=command,binary_sha256=sha(b),
            initial_field_sha256=sha(native/'field.csv'),OMP_NUM_THREADS=threads,OMP_DYNAMIC=False))
        start=perf_counter()
        with (dest/'stdout.log').open('w') as stream:
            subprocess.run(command,cwd=dest,env=environment,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1800)
        a=analyze(dest)
        cases[label]=dict(analysis=a,elapsed_s=perf_counter()-start,
            checks=native_checks(a,carrier['native_target_velocity_m_s']))
        print(label,cases[label]['elapsed_s'],cases[label]['checks'],flush=True)
    fields={label:np.genfromtxt(base/label/'field.csv',delimiter=',',names=True) for label in cases}
    if not np.array_equal(fields['serial'][['i','j','k']],fields['omp'][['i','j','k']]):
        raise ValueError('Parallel geometry/indexing differs.')
    difference=max(float(np.max(np.abs(fields['serial'][key]-fields['omp'][key]))) for key in ('ux_m_s','uy_m_s','uz_m_s'))
    peak=max(float(np.max(np.abs(fields['serial'][key]))) for key in ('ux_m_s','uy_m_s','uz_m_s'))
    checks=dict(all_native=all(all(c['checks'].values()) for c in cases.values()),
        pointwise_velocity=difference/peak<=1e-8,
        superficial_velocity=abs(cases['omp']['analysis']['superficial_velocity_m_s']/cases['serial']['analysis']['superficial_velocity_m_s']-1)<=1e-10)
    result=dict(schema='pore2process.stock-openmp-crosscheck.v1',run_id=rid,source_sha256=sources,
        hydraulic_summary=str((hydraulic/'summary.json').relative_to(ROOT)),
        hydraulic_summary_sha256=sha(hydraulic/'summary.json'),
        cases=cases,threads=args.threads,build=build,binary=str(binary.resolve()),binary_sha256=sha(binary),
        max_absolute_velocity_difference_m_s=difference,max_relative_velocity_difference=difference/peak,
        measured_speedup=cases['serial']['elapsed_s']/cases['omp']['elapsed_s'],
        checks=checks,checks_passed=all(checks.values()),vendor=verify_vendor(),
        native_runs=2,new_DWSIM_runs=0,new_adsorption_predictions=0)
    result['artifacts_sha256']={str(p.relative_to(base)):sha(p) for p in base.rglob('*') if p.is_file()}
    write_json(base/'summary.json',result)
    write_json(OUTPUT/'omp_latest.json',dict(summary=str((base/'summary.json').relative_to(ROOT))))
    print('OpenMP cross-check result',checks,'speedup',result['measured_speedup'],flush=True)
    if not result['checks_passed']:raise ValueError('Stock OpenMP native cross-check failed.')


if __name__=='__main__':main()
