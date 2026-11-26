"""Uniform salinity under transverse gravity should establish a pressure head."""
from pathlib import Path
import os
import subprocess

import numpy as np

from pore2process_demo.io import read_json, write_json
from .run import OUT


def main():
    source=OUT/"N8_low_down"
    c=read_json(source/"inputs.json")
    args=(source/"native_args.txt").read_text().splitlines()
    args[3]="1"
    out=OUT/"uniform_salt"
    out.mkdir(exist_ok=True)
    result=[]
    with (out/"native.log").open("w") as log:
        with subprocess.Popen(args,cwd=out,stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,text=True,bufsize=1,
                              env={**os.environ,"OMP_NUM_THREADS":"2","OMP_DYNAMIC":"FALSE"}) as proc:
            def wait_for(tag):
                for line in proc.stdout:
                    log.write(line); log.flush()
                    if line.startswith(tag+" "):
                        return list(map(float,line.split()[1:]))
                raise RuntimeError("Native hydrostatic control terminated early")

            def solve(fraction,steps):
                np.full(c["native_rows"],fraction,dtype="<f8").tofile(out/"salt.bin")
                proc.stdin.write(f'{steps} {c["gravity_accel"]} 1\n'); proc.stdin.flush()
                raw=wait_for("SALT_DONE")
                return raw,np.fromfile(out/"water.bin",dtype="<f8").reshape(-1,3)

            wait_for("SALT_READY")
            raw,base=solve(0.,2500)
            for f in (0.1,1.0):
                old=base
                for sweeps in (600,1200,2400,4800):
                    # The increments double after the first two samples.
                    raw,u=solve(f,sweeps if sweeps==600 else sweeps//2)
                    scale=float(np.linalg.norm(base))
                    result.append(dict(feed_fraction=f,total_sweeps=sweeps,
                        mean_z_m_s=float(u[:,2].mean()),baseline_mean_z_m_s=float(base[:,2].mean()),
                        fraction_mean_axial_speed=float((u[:,2]-base[:,2]).mean()/base[:,0].mean()),
                        rms_velocity_change_relative=float(np.linalg.norm(u-base)/scale),
                        step_rms_change_relative=float(np.linalg.norm(u-old)/scale),
                        native_Q_m3_s=raw[0],max_Mach=raw[1],max_density_departure=raw[2]))
                    old=u
                    if sweeps>=1200 and result[-1]["step_rms_change_relative"]<1e-5:
                        break
            proc.stdin.write("-1 0 1\n"); proc.stdin.flush()
            if proc.wait()!=0:
                raise RuntimeError("Native hydrostatic control failed")
    write_json(out/"control.json",dict(rows=result,baseline_mean_axial_m_s=float(base[:,0].mean()),
               baseline_water_sha256=c["water_sha256"],units="velocity m/s, discharge m3/s"))


if __name__=="__main__":
    main()
