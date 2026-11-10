"""Render accepted numerical findings into the existing methods report."""
from pathlib import Path
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read


def update(s,path):
    path=Path(path);n=read(path.parent/'native_summary.json');v=read(path.parent/'verification.json')
    if s['status']!='passed' or v['status']!='passed':raise ValueError('Unverified results cannot update the report.')
    rows=s['comparison'];nom=[r for r in rows if r['inlet_Re']==1388 and r['flux_LMH']==12 and r['fixed_multiplier']==1 and r['viscosity_multiplier']==1]
    lines=['## Verified results','',f"Accepted native run: `{s['native_path']}`.",f"Accepted process run: `{path}`.",'',
        f"The suite passes **{len(n['cases'])} native D3Q27/TRT cases** and **{len(s['cases'])+1} actual DWSIM solve/save/reload cases**, including the independent seed check. A one-iteration recycle is separately rejected.",'',
        '| Native quantity | Result |','| --- | ---: |',
        f"| Accepted C | {n['geometry_coefficient_m4']:.12e} m4 |",
        f"| Measured Re envelope | {n['verified_Re_range'][0]:.6f}–{n['verified_Re_range'][1]:.6f} |",
        f"| 64/48 spatial-pair C change | {100*(n['cases']['grid64']['geometry_coefficient_m4']/n['cases']['grid48']['geometry_coefficient_m4']-1):.4f}% |",
        f"| Finest C versus independent Julia reference | {100*(n['cases']['grid64']['geometry_coefficient_m4']/n['reference_C_m4']-1):.4f}% |",
        f"| Maximum lattice Mach | {max(c['raw']['max_Mach'] for c in n['cases'].values()):.6f} |",
        f"| Maximum density departure | {max(c['raw']['max_density_departure'] for c in n['cases'].values()):.6e} |",
        f"| Time-step C change | {100*abs(n['cases']['time']['geometry_coefficient_m4']/n['cases']['high1388_n64']['geometry_coefficient_m4']-1):.6f}% |",
        f"| Domain C change | {100*abs(n['cases']['domain']['geometry_coefficient_m4']/n['cases']['high1388_n64']['geometry_coefficient_m4']-1):.6f}% |",'',
        'The prototype protocol contains a nominal tau; actual saved case requests/raw files are authoritative. Baselines use tau=0.8, operating points use tau=0.5+N/64000, and the temporal check uses tau=0.5005. All use magic=3/16. Neither a source-paper operating regime nor physical transition is inferred from a stationary periodic laminar branch.','',
        '### Material comparison at fixed 12 LMH and inlet Re=1388','',
        'Recovery is 0.5, pump efficiencies 0.75, temperature 25 °C and fixed-impedance multiplier 1.','',
        '| Support interior | Mean TMP (kPa) | Shaft energy (kWh/m3) | Saving versus 4 µm | Fixed resistance share | Circulation power share |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for r in nom:
        label={'straight':'4 µm','widened_5um':'5 µm','widened_6um':'6 µm'}[r['support']]
        lines.append(f"| {label} | {r['mean_TMP_Pa']/1000:.6f} | {r['specific_shaft_energy_kWh_m3']:.9f} | {100*r['support_energy_saving_fraction']:.5f}% | {100*r['fixed_resistance_fraction']:.4f}% | {100*r['circulation_power_fraction']:.4f}% |")
    lines += ['', 'Support geometry still lowers the transmembrane load. At these circulation rates, most pump power supplies the axial channel loss; transmembrane resistance dominance and total process-energy dominance are different questions. The historical low-Re example retains its own descriptor/discretization evidence and is not treated as a strict controlled comparison with this new study.','',
        '| Support interior | Fixed-impedance multiplier at equal support/fixed contributions |','| --- | ---: |']
    for r in nom:
        c=next(c for c in s['cases'] if c['id']==r['case_id']);a=c['translation']['axial']
        lines.append(f"| {r['support']} | {a['support_impedance_Pa_s_m']/a['fixed_impedance_Pa_s_m']:.6f} |")
    lines += ['', 'Above those multipliers the fixed contribution dominates; below them the support dominates in this series law. The sampled 0.5/1/2 multipliers are above all three crossings. These are effective non-support impedances, not separately measured intrinsic selective-layer resistances.','',
        '### Conditional operating selections','',
        '| Support | Mean flux (LMH) | Fixed factor | Mapping | Capacity fraction | Policy | Selected Re | Shaft energy (kWh/m3) |',
        '| --- | ---: | ---: | --- | ---: | --- | ---: | ---: |']
    examples=[('straight',12,1,'equal_Re',.95,'linear_nominal'),('straight',12,1,'equal_velocity',.95,'linear_nominal'),
              ('widened_6um',15,.5,'equal_Re',.8,'linear_nominal'),('widened_6um',15,.5,'equal_Re',.8,'monotonic_lower')]
    for sup,j,z,m,cap,policy in examples:
        q=next(q for q in s['selections'] if (q['support'],q['mean_flux_LMH'],q['fixed_multiplier'],q['mapping'],q['capacity_fraction'],q['policy'],q['transfer_factor'])==(sup,j,z,m,cap,policy,1.))
        rr='No certified grid point' if q['selected_Re'] is None else f"{q['selected_Re']:g}"
        energy='Unavailable' if q['specific_shaft_energy_kWh_m3'] is None else f"{q['specific_shaft_energy_kWh_m3']:.9f}"
        lines.append(f'| {sup} | {j} | {z} | {m} | {cap} | {policy} | {rr} | {energy} |')
    lines += ['', 'These are conditional discrete-grid selections. A missing certified point can reflect unavailable source coverage or indeterminate monotonic bounds; it is not proof of physical infeasibility. Transfer factor is 1 in this table. Full policies and transfer sensitivities are in `selections.csv`.','',
        '### Hydraulic-property and quality sensitivity','',
        '| Inlet Re | 1.2× viscosity energy / water-anchor energy, 6 µm support at 12 LMH |','| --- | ---: |']
    for re in s['protocol']['process']['Re_targets']:
        rr=[r for r in rows if r['support']=='widened_6um' and r['flux_LMH']==12 and r['fixed_multiplier']==1 and r['inlet_Re']==re]
        low=next(r for r in rr if r['viscosity_multiplier']==1);high=next(r for r in rr if r['viscosity_multiplier']==1.2)
        lines.append(f"| {re} | {high['specific_shaft_energy_kWh_m3']/low['specific_shaft_energy_kWh_m3']:.6f} |")
    lines += ['', 'This sensitivity keeps Re fixed: increasing assumed hydraulic viscosity also increases channel velocity/flow. Water-carrier DWSIM thermodynamics are retained; this is not measured biological-feed rheology.','']
    c=next(c for c in s['cases'] if c['id']==nom[0]['case_id']);t=c['translation']
    quality=[(q['assumed_local_rejection'],q['balance']['permeate_COD_kg_m3']*1000) for q in t['quality_sensitivity']]
    lines += [f"At the nominal point, illustrative rejection 0.65 gives permeate COD {t['feed']['permeate_COD_kg_m3']*1000:.6f} mg/L. Chosen retention sensitivities give {quality[0][1]:.6f} mg/L at rejection {quality[0][0]} and {quality[1][1]:.6f} mg/L at rejection {quality[1][0]}. These are external conservative balances on the same hydraulics; no new DWSIM species or reaction runs are counted.",'',
        '| Omitted source point Re | Leave-one-interior-point-out relative error |','| --- | ---: |']
    for d in s['source_shape_diagnostics']:lines.append(f"| {d['held_out_source_Re']:g} | {100*d['relative_error']:.4f}% |")
    lines += ['', 'These errors diagnose interpolation shape within the same four-point series. They are not independent validation or experimental uncertainty.','',
        '### Process and reduced-model checks','',
        f"Maximum unit/overall relative mass residual: {max(c['validation']['max_mass_relative_residual'] for c in s['cases']+[s['seed_validation']]):.6e}.",
        f"Maximum enthalpy residual normalized by total shaft power: {max(c['validation']['max_enthalpy_relative_residual'] for c in s['cases']+[s['seed_validation']]):.6e}.",
        f"Largest omitted-axial-momentum/friction upper screen: {100*max(c['translation']['axial']['neglected_axial_momentum_to_friction_upper_screen'] for c in s['cases']):.6f}% (limit 2%).",'',
        f"The checker replays {v['capacity_decisions']} source-capacity comparisons and {v['discrete_policy_selections']} discrete policy selections. Independent upstream replay recorded: `{v['upstream_rechecked']}`.",'',
        '[Scientific figure](operating_map.png) · [Compact numerical reference](operating_map_reference.json) · [Reproduction commands](../operating_map/README.md)','',
        'The standalone map adds no new Panel browser controls. The existing low-Re interface retains its previous browser verification; no new browser pass is claimed for this CLI/report extension. Reacting and non-isothermal additions remain deferred.','']
    file=ROOT/'docs/operating_window.md';text=file.read_text().split('## Verified results')[0]
    start='Implementation record, 9 November 2026. **Verification is in progress.**\nFinal accepted run identities and results will replace this status after the\nnative suite and actual DWSIM matrix complete.'
    text=text.replace(start,'Implementation and verification record, 9–10 November 2026. The bounded native\nand actual-process checks pass; scientific transfer limits remain explicit below.')
    text=text.replace('Acceptance results, figure and run identities will be added after verification.','Accepted results and run identities follow below.')
    file.write_text(text.rstrip()+'\n\n'+'\n'.join(lines))
