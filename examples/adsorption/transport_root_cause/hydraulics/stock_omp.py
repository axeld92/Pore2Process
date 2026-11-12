"""External build of an upstream OpenMP CPU_SISD case and matching core.

This helper does not change the vendor tree or global OpenMP environment.
The caller supplies project application, build and record directories and
sets OMP_NUM_THREADS explicitly for each native process.
"""
from pathlib import Path
import shutil
import subprocess
from time import perf_counter

from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from pore2process_demo.runtime.openlb import OPENLB,verify_vendor


def build_case_omp(source_directory,build_directory,output_directory):
    source,build,out=map(lambda p:Path(p).resolve(),
        (source_directory,build_directory,output_directory))
    for path in (source,build,out):
        if path==OPENLB.resolve() or OPENLB.resolve() in path.parents:
            raise ValueError('Application, build and output must remain outside the vendor tree.')
    vendor=verify_vendor()
    build.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    core=build/'core';core.mkdir(exist_ok=True)
    source_sha={}
    for path in source.iterdir():
        if path.is_file():
            source_sha[path.name]=sha(path);shutil.copy2(path,build/path.name)
    config=build/'config.mk'
    text=config.read_text()
    if 'PARALLEL_MODE := OFF' in text:
        text=text.replace('PARALLEL_MODE := OFF','PARALLEL_MODE := OMP')
    elif 'PARALLEL_MODE := OMP' not in text:
        raise ValueError('Expected explicit OFF or OMP parallel configuration.')
    config.write_text(text+'\nOMPFLAGS := -fopenmp\n')
    compiler='/opt/homebrew/bin/g++-15'
    flags=['-O2','-Wall','-std=c++20','-pthread','-fPIC','-DPLATFORM_CPU_SISD',
        '-DPARALLEL_MODE_OMP','-fopenmp','-DOLB_VERSION="1.9r0"',
        f'-I{OPENLB/"src"}','-I/opt/homebrew/include']
    native_sources=('communication/mpiManager.cpp','communication/ompManager.cpp',
        'core/olbInit.cpp','core/expr.cpp','io/ostreamManager.cpp')
    commands=[];objects=[];start=perf_counter()
    with (out/'build.log').open('w') as log:
        for relative in native_sources:
            obj=core/(relative.replace('/','_')+'.o');objects.append(obj)
            command=[compiler,*flags,'-c',str(OPENLB/'src'/relative),'-o',str(obj)]
            commands.append(command);subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        command=['ar','rcs',str(core/'libolbcore.a'),*map(str,objects)]
        commands.append(command);subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
        (build/'capillary3d.o').unlink(missing_ok=True)
        command=['make','-j2','GIT_VERSION=1.9r0',f'OLB_ROOT={OPENLB}',f'PROJECT_ROOT={ROOT}','onlysample']
        commands.append(command);subprocess.run(command,cwd=build,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    metadata=dict(commands=commands,source_sha256=source_sha,core_flags=flags,
        core_sources=list(native_sources),external_config=config.read_text(),
        build_s=perf_counter()-start,compiler=subprocess.check_output([compiler,'--version'],text=True),
        vendor_integrity=vendor,binary=str(build/'capillary3d'),binary_sha256=sha(build/'capillary3d'))
    verify_vendor();write_json(out/'build.json',metadata)
    return metadata
