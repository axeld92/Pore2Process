# Equilibrium sensitivity at the Figure 7 particle sizes

Test whether the available equilibrium evidence supports transfer between
0.80 and 0.45 mm particles at 43.8 mL/min, 150 mg/L phenol and a 10 cm bed.
Use the saved N160 TRT tracer models at the two established section Péclet
numbers. Keep the bed liquid volume and the density-derived 7.032455 g
dry carbon mass fixed for both sizes. Reuse the original per-size
film correlation, contact shares, travel-time memory and mixed joins.
The tracer geometry and periodic section mapping are synthetic.

Use the prior half-direct, half-slow carbon scenario with a common
added resistance. Derive it from the archived 0.80 mm film and 1%
slow/direct conductance ratio. Apply
`k_slow/k_film = 1/(1 + R k_film)` at both particle sizes, as in the
earlier bed-depth transfer. Retain the full dry mass. Film, access and
transport inputs stay fixed across the equilibrium cases.
[Travel-time qualification](../disordered_travel_memory/REPORT.md)
· [Access resistance and scale-up](../carbon_access_transfer/REPORT.md)

Use four predeclared material cases.

- Alhamed Table 1 Langmuir with the archived stirred-bath kinetics.
- The recipe-matched GAC4 103 mg/g Langmuir plateau with the previously
  checked affinity and the same Langmuir kinetics.
- The 50/50 loading blend chosen on the 10 and 15 cm Figure 5 cases at
  23.3 mL/min, with its batch-qualified diffusivity.
- Alhamed Table 1 Freundlich with its batch-qualified diffusivity.

The GAC4 plateau is an extrapolated fit and the affinity comes from the
2009 batch endpoint consistency check. The four cases are sensitivities
for the same material and inventory. Their range is conditional on
the tested equations and kinetics.
[Material checks](../breakthrough_hypotheses/batch/REPORT.md)
· [Frozen blend selection](../breakthrough_hypotheses/batch/conditional/README.md)

Freeze all predictions and refinement results before scoring the two
Figure 7 curves. Check conservation, physical inventory, interface
flux and time/radial refinement. Use the existing 25% 5%-crossing,
20% midpoint and outlet/feed RMSE 0.10 limits. Compare retained phenol
at the last marker in each series with its digitization reading range.
Report post-midpoint retention using that series' published midpoint.
No Figure 7 point determines a model input.

Use the saved, verified DWSIM pump powers at the same flow and dry carbon
mass for the process comparison. The archived synthetic FCC hydraulic
bed is 14.325 cm long at 0.60 porosity. The experimental adsorption bed is
10 cm long at the reported 0.427 porosity. With the 5% threshold, the usable
volume is `Q t_05` and the bed-friction shaft energy is `P_shaft t_05`.
Report these quantities for the measured curve crossings, including
the existing digitization sensitivity, and for each simulated case.
The measured-crossing calculation combines experimental service time
with simulated pressure loss. The same material and transport assumptions
must pass both sizes for the predicted service-time comparison.
[Alhamed Figure 7](https://doi.org/10.1016/j.jhazmat.2009.05.002)
· [DWSIM pump results](../column_transport/REPORT.md)
