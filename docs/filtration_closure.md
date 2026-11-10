# Filtration hydraulic milestone closure

Status: closed for now on 10 November 2026. The existing capillary-network,
low-Re filtration interface and prior connected examples remain available.
The final higher-Re operating map completes this hydraulic example; isothermal
adsorption is the next case in development and is published separately.

## Engineering conclusion

Improving a porous component does not necessarily reduce the whole process's
energy by the same fraction. At fixed 12 LMH, inlet Re=1388, 50% recovery and
unchanged effective non-support impedance, widening the support interior from
4 to 6 µm lowers mean TMP from 1.2211 to 1.0831 kPa. Total pump shaft energy
falls from 0.056703822 to 0.056600813 kWh/m³: a 0.18166% saving. Circulation
accounts for 98.7234% of baseline shaft duty.

The fixed effective contribution rises from 68.5026% to 77.2330% of membrane
resistance. Those percentages describe the membrane resistance partition;
they are distinct from the allocation of process pump duty. The comparison
holds filtration flux, circulation, recovery and the fixed impedance constant.

## Completed evidence

- 16 native D3Q27/TRT cases, with 144 numerical screens and a measured
  fully developed laminar-branch envelope near Re=1000–2000.
- 91 actual DWSIM solve/save/reload cases, including a perturbed-seed check;
  an insufficient-iteration case is separately rejected.
- 902 Julia assertions and 74 filtration-only Python tests recorded at closure.
- Replay of 682 local artifacts, 1,620 capacity comparisons and 1,458 discrete
  policy selections, with the original simulator source audit retained.

The measured native coefficient drives Julia's distributed model and actual
DWSIM port pressures and pump duties. Native hydrodynamic verification and
conditional literature-capacity transfer are documented separately. The
higher-Re map adds no new Panel controls or browser-verification claim.

## Preserved scope and next work

The current map evaluates declared interpolation, monotonic-bound and
Reynolds/velocity-transfer policies. It does not establish physical transition
or a measured fouling-time law. The scientific report retains all source,
scenario, operating-range and reduced-model assumptions.

Further filtration physics is optional future work. The next active research
question concerns particle size, treatment capacity and pumping energy in an
isothermal adsorption bed. This publication includes only the closed filtration
source and evidence; concurrent adsorption files and detailed claims remain
outside the snapshot.

[Final scientific report](operating_window.md) ·
[Ordered reproduction](reproduction.md) ·
[Operating-map commands](../examples/filtration/operating_map/README.md)
