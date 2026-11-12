"""Read-only checks of the actual initial scalar phase in archived populations."""
import json
from pathlib import Path
import numpy as np
from pore2process_demo.paths import ROOT
from pore2process_demo.io import sha,write_json
from .run import OUTPUT,HERE,check


def initial_profile(directory,phase_x_m):
    directory=Path(directory);raw=json.loads((directory/'raw.json').read_text())
    a=np.genfromtxt(directory/'audit_1.csv',delimiter=',',names=True)
    x=(a['i']+.5)*raw['dx_m']
    exact=.5+.25*np.sin(2*np.pi*(x-phase_x_m)/raw['side_m'])
    error=float(np.max(abs(a['concentration_before']-exact)))
    fluid=a['material']!=2;initial=a['concentration_before'][fluid]
    return dict(phase_x_m=phase_x_m,nodes_checked=len(a),maximum_initial_profile_error_kg_m3=error,
                initial_liquid_mean_C_kg_m3=float(np.mean(initial)),initial_liquid_SD_C_kg_m3=float(np.std(initial)),
                initial_liquid_inventory_kg=float(np.sum(initial)*raw['dx_m']**3),
                physical_liquid_volume_m3=float(np.sum(fluid)*raw['dx_m']**3),
                tolerance_kg_m3=1e-12,passed=error<=1e-12)


def main():
    verification=check();summary=ROOT/verification['summary'];s=json.loads(summary.read_text())
    baseline=ROOT/s['baseline_summary']
    checks=dict(unshifted=initial_profile(baseline.parent/'native/n32_sine',0.),
        shifted=initial_profile(summary.parent/'native/shifted_sine',s['phase_x_m']))
    if not all(c['passed'] for c in checks.values()):raise ValueError('Actual initialized scalar phase differs.')
    normalized={'unshifted':s['unshifted']['result']['periodic_wave_final']['sd_C_kg_m3']/checks['unshifted']['initial_liquid_SD_C_kg_m3']}
    normalized['shifted']=s['shifted']['shifted_sine']['result']['periodic_wave_final']['sd_C_kg_m3']/checks['shifted']['initial_liquid_SD_C_kg_m3']
    result=dict(summary=str(summary.relative_to(ROOT)),summary_sha256=sha(summary),
                checker_sha256=sha(Path(__file__)),checks=checks,final_SD_over_initial_SD=normalized,
                relative_change_in_normalized_SD=normalized['shifted']/normalized['unshifted']-1,
                new_native_runs=0,new_DWSIM_runs=0)
    write_json(OUTPUT/'initial_profile_checks.json',result)
    write_json(HERE/'initial_profile_reference.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
