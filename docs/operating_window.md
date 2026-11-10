# Higher-Re verification and conditional filtration operating map

Implementation and verification record, 9–10 November 2026. The bounded native
and actual-process checks pass; scientific transfer limits remain explicit below.

## Research question and boundaries

At fixed mean filtration flux, unchanged fixed membrane impedance and specified
recovery, how do support geometry and circulation change total pumping energy?
Which design/operation pairs meet a declared literature-informed capacity
screen, and how sensitive is that decision to the transfer assumptions?

The model is isothermal clean initial hydraulics. The synthetic channel is
horizontal, gravity is neglected, and external permeate pressure is uniform.
The paper's rig is vertically mounted; no source hydrostatic profile is inferred. A short-step critical flux
is a protocol-dependent observation, not an indefinitely sustainable operating
flux. Physical fouling time, measured effluent viscosity and same-condition
retention versus circulation remain unavailable. Reaction and non-isothermal
extensions are deferred.

## Native implementation and acceptance

All project application source is in `operating_map/openlb/`; external build
products are in `.local/operating_map/build_trt27/`. Original OpenLB files are
not modified. The installed 1.9r0 examples and stock components determine the
API and lifecycle: mesh, geometry, lattice, initial conditions, simulation,
results. See the official [example structure](https://openlb.gitlab.io/release/content/11_step_by_step.html)
and [forcing API](https://openlb.gitlab.io/release/api/d2/d86/namespaceolb_1_1forcing.html).

The application uses the existing D3Q27 descriptor, TRT collision,
Kupershtokh forcing, halfway bounce-back, STL/SDF geometry, plane-integral
functors and ValueTracer. The alias composes upstream types; no collision,
velocity set, geometry class, data structure or numerical library is recreated.
TRT magic is held at 3/16 for all checks. Warm starts use actual native velocity and non-equilibrium stress fields,
never Julia's independent reference as an initializer. OpenLB's stock
`iniRegularized` initializes the latter after lattice-unit conversion; older
velocity-only source files use stock equilibrium initialization. Full-profile
initial-estimate invariance is checked separately.

Both discharge and mean lattice kinetic energy must become stationary under
stock six-sample tracers (relative criterion 1e-8). Cases also observe at least
0.05 hydraulic-diameter viscous time units. Checks retain grid, time-step,
periodic-domain, initial-estimate, Mach, density and mass-conservation screens.
Reported Reynolds numbers come from the measured plane flux, declared area
and hydraulic diameter; the converter's reference velocity is not the imposed
channel velocity. Raw non-equilibrium moment columns are internal warm-start
data, not an independent physical wall-stress measurement.
The circular channel is checked against Hagen–Poiseuille; the synthetic star
against an independent Julia polygon Poisson reference. Pressure scaling is
checked on the higher-Re native cases before the coefficient is used in Julia.

The first execution ceiling was 3,600 seconds. The longer-domain case was
restarted from its own finite native velocity/stress field with a 7,200-second
allowance. Its earlier attempt and field provenance are retained. This changes
the initial estimate and execution budget, with the same physical case,
minimum observation interval and 1e-8 stationarity criteria.

These are fully developed **laminar-branch** checks in a straight periodic
window. They do not verify entry losses, physical transition, turbulence or a
full-length permeating-wall CFD model. The 0.3 m distributed model remains a
slowly varying approximation with its earlier inertia/suction screens intact.
An additional maximum-local-flux screen bounds the omitted axial momentum
term relative to friction by `2*beta*rho*C*P*Jmax/(mu*A^2)`, with
`beta <= umax/umean` obtained from the measured native field. It must remain
below 2%. This is an a priori reduced-model screen, not a permeating-wall CFD
validation.
Mean wall traction derived from `G*A/P` is an integral momentum-balance value,
not a resolved local wall-shear distribution.

## Published evidence and conditional transfer

The auditable extraction is [operating_window_evidence.json](../literature/operating_window_evidence.json).
The original four-paper extraction and its source hashes are retained.

- [Chiu et al. (2006)](https://doi.org/10.1016/j.seppur.2005.11.028), §3.2/Fig. 6:
  biological-process effluent, nominal 0.2 µm membrane, pH 7.5, 20–25 °C.
  Source Reynolds points are 1388, 2776, 3701 and 4626. Written flux anchors
  are 21.4 and 72.5 LMH; the other fluxes are digitized with reading bounds.
  The paper calls Re=1388 laminar and Re=2776 turbulent. The curve is not a
  laminar constitutive law. Its 30/15-minute TMP-step criterion is preserved.
- [Chiu and James (2005)](https://doi.org/10.1016/j.memsci.2005.01.037) and
  [Garcia and Chiu (2008)](https://doi.org/10.1016/j.memsci.2008.08.033) concern
  TiO2 suspensions. Their critical-flux data overlap; they are not independent
  biological-feed validation. Garcia's absolute bundle-flow convention remains
  unresolved and is not used to define the synthetic channel flow.
- [Frederic et al. (2018)](https://doi.org/10.1016/j.memsci.2018.09.021)
  supplies a porous cross-section model. Its arbitrary fouling coefficient
  does not supply measured crossflow-dependent capacity or a physical time fit.

The fixed hydraulic impedance is transferred from the earlier **paired Zhu
specimen** extraction; Chiu uses another material/feed. Combining them is a
scenario, not a matched experimental prediction. No layer resistance is inferred
by subtracting unrelated specimens.

Julia evaluates an explicit piecewise-linear interpolation within the source
range and separate shape-free monotonic brackets. There is no extrapolation.
A leave-one-interior-point-out calculation diagnoses interpolation shape; it
is not independent validation. The ±2-pixel reading bounds are sensitivity choices, not experimental error
bars. The authors' reported maximum 5% reproducibility is retained separately;
it is not a confidence interval or used to invent point uncertainties.

Equal-Re and equal-velocity transfer are separate policies. Equal velocity uses
only the source's reported pair (Re=1388, u=0.6 m/s); it does not infer measured
feed viscosity. Capacity fractions 0.8/0.95/1 and transfer factors 0.8/1/1.2 are
chosen scenarios. A conservative local screen compares **maximum local flux**
with the transferred source capacity. Using a module-average experimental Jcrit
locally is itself a declared assumption. Monotonic bounds may leave a decision
indeterminate even when the linear scenario permits it.

COD rejection remains the illustrative nominal constant 0.65, with external
0.5/0.8 balance sensitivities on each same hydraulic result; no circulation-dependent
quality improvement is invented. The external tracer balance records module,
fresh, permeate and purge concentrations. COD removal means diversion to purge,
not reaction. The carrier in DWSIM is water. These quality scenarios do not
constitute additional DWSIM species simulations or a measured quality-benefit
law. No regulatory permeate-quality threshold is invented for the selection.

## Design map and actual process calculation

The evaluated grid is three supports (4/5/6 µm interior; fixed 4 µm openings),
mean fluxes 8/12/15 LMH, inlet Re=1388/1600/1850 and fixed-impedance multipliers
0.5/1/2. Within each geometry comparison, the flux, circulation, recovery,
viscosity and fixed impedance stay identical. Nine additional hydraulic-viscosity-multiplier 1.2 scenarios
at 12 LMH/fixed factor 1 retain the specified inlet Re; therefore channel speed
and flow rise with viscosity. They keep water-carrier thermodynamics in DWSIM
and represent an assumed hydraulic-property sensitivity, not measured effluent
rheology. Main policy selection uses only multiplier 1. Unphysical cases are
rejected explicitly rather than repaired.

Julia receives the accepted measured OpenLB coefficient and existing native
TRT support resistances. `reference_membrane` retains the original calibration
anchor at its reference flux; the current case flows and pressures are in
`axial` and `process`. Module length (0.3 m) and module-inlet COD are explicitly
retained in the process input snapshot. The existing distributed channel calculation supplies
actual port loads to the stock DWSIM recycle. Two pump shaft powers, cooling
heat, actual mass/enthalpy balances, phase guards, tear convergence, saved
connectivity and reload invariance are checked. The stock forced-liquid IAPWS
route and the previous strict acceptance tolerances are retained. Only the
iteration ceiling is raised for the high circulation ratio; the nominal tear
mass seed uses the exact overall mass balance, with a +2 K initial temperature
perturbation to require a genuine stock recycle iteration. A mass-only exact
seed stopped at zero iterations with a cooling-duty mismatch and was rejected. A half-flow seed at a -1 K initial temperature offset, with tighter stop
criteria, separately checks invariance, and a one-iteration recycle must fail.

The result is pump shaft energy per reference filtrate volume. It excludes
motor/cooling electricity, cleaning, capital and long-term production downtime.
Selections are minima **on the evaluated discrete grid under stated policies**;
they are not validated physical or continuous optima.

## Reproduction and verification record

See [exact commands and resumption](../operating_map/README.md). Existing locks
and runtimes are reused; no dependency installation or system changes are
required. `results/operating_map/preflight.json` records the actual machine.

The direct Julia suite passes **902 assertions**, and the full Python suite
passes **74 tests** in the normal host environment. The latter replays the
real native/Julia/DWSIM reference and checks source bounds, unavailable values,
invalid requests and changed artifacts. Saved-record tests are distinct from
the 91 fresh external solves performed by the study. The scientific figure
was visually inspected. Existing dependency locks are unchanged.

The first process requests exposed our adapter's old 5,000-iteration input
guard; it was extended to 20,000 without changing stock DWSIM or convergence
tolerances, and a fresh source-snapshotted run passed. Its independent
half-flow/-1 K seed check takes 8,796 iterations and agrees with nominal duties
within 1e-8 relative. Earlier failures remain saved separately.

Accepted results and run identities follow below.
The original low-Re Panel showcase and completed examples remain available.

## Verified results

Accepted native run: `/Users/dula/Documents/ChatGPT/Pore2Process/results/operating_map/native/runs/20261109T205248Z-d6fa71dd/summary.json`.
Accepted process run: `/Users/dula/Documents/ChatGPT/Pore2Process/results/operating_map/process/runs/20261109T233235Z-1baf9ec0/summary.json`.

The suite passes **16 native D3Q27/TRT cases** and **91 actual DWSIM solve/save/reload cases**, including the independent seed check. A one-iteration recycle is separately rejected.

| Native quantity | Result |
| --- | ---: |
| Accepted C | 2.626622470254e-12 m4 |
| Measured Re envelope | 999.999697–1999.999393 |
| 64/48 spatial-pair C change | 1.4358% |
| Finest C versus independent Julia reference | -1.4796% |
| Maximum lattice Mach | 0.077399 |
| Maximum density departure | 2.775558e-15 |
| Time-step C change | 0.000031% |
| Domain C change | 0.000000% |

The prototype protocol contains a nominal tau; actual saved case requests/raw files are authoritative. Baselines use tau=0.8, operating points use tau=0.5+N/64000, and the temporal check uses tau=0.5005. All use magic=3/16. Neither a source-paper operating regime nor physical transition is inferred from a stationary periodic laminar branch.

### Material comparison at fixed 12 LMH and inlet Re=1388

Recovery is 0.5, pump efficiencies 0.75, temperature 25 °C and fixed-impedance multiplier 1.

| Support interior | Mean TMP (kPa) | Shaft energy (kWh/m3) | Saving versus 4 µm | Fixed resistance share | Circulation power share |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4 µm | 1.221104 | 0.056703822 | 0.00000% | 68.5026% | 98.7234% |
| 5 µm | 1.123836 | 0.056631254 | 0.12798% | 74.4316% | 98.8489% |
| 6 µm | 1.083071 | 0.056600813 | 0.18166% | 77.2330% | 98.9017% |

Support geometry still lowers the transmembrane load. At these circulation rates, most pump power supplies the axial channel loss; transmembrane resistance dominance and total process-energy dominance are different questions. The historical low-Re example retains its own descriptor/discretization evidence and is not treated as a strict controlled comparison with this new study.

| Support interior | Fixed-impedance multiplier at equal support/fixed contributions |
| --- | ---: |
| straight | 0.459798 |
| widened_5um | 0.343516 |
| widened_6um | 0.294783 |

Above those multipliers the fixed contribution dominates; below them the support dominates in this series law. The sampled 0.5/1/2 multipliers are above all three crossings. These are effective non-support impedances, not separately measured intrinsic selective-layer resistances.

### Conditional operating selections

| Support | Mean flux (LMH) | Fixed factor | Mapping | Capacity fraction | Policy | Selected Re | Shaft energy (kWh/m3) |
| --- | ---: | ---: | --- | ---: | --- | ---: | ---: |
| straight | 12 | 1 | equal_Re | 0.95 | linear_nominal | 1388 | 0.056703822 |
| straight | 12 | 1 | equal_velocity | 0.95 | linear_nominal | 1850 | 0.100151975 |
| widened_6um | 15 | 0.5 | equal_Re | 0.8 | linear_nominal | 1600 | 0.059903053 |
| widened_6um | 15 | 0.5 | equal_Re | 0.8 | monotonic_lower | No certified grid point | Unavailable |

These are conditional discrete-grid selections. A missing certified point can reflect unavailable source coverage or indeterminate monotonic bounds; it is not proof of physical infeasibility. Transfer factor is 1 in this table. Full policies and transfer sensitivities are in `selections.csv`.

### Hydraulic-property and quality sensitivity

| Inlet Re | 1.2× viscosity energy / water-anchor energy, 6 µm support at 12 LMH |
| --- | ---: |
| 1388 | 1.721894 |
| 1600 | 1.723551 |
| 1850 | 1.724806 |

This sensitivity keeps Re fixed: increasing assumed hydraulic viscosity also increases channel velocity/flow. Water-carrier DWSIM thermodynamics are retained; this is not measured biological-feed rheology.

At the nominal point, illustrative rejection 0.65 gives permeate COD 24.070272 mg/L. Chosen retention sensitivities give 34.377767 mg/L at rejection 0.5 and 13.757776 mg/L at rejection 0.8. These are external conservative balances on the same hydraulics; no new DWSIM species or reaction runs are counted.

| Omitted source point Re | Leave-one-interior-point-out relative error |
| --- | ---: |
| 2776 | -9.8635% |
| 3701 | -3.1495% |

These errors diagnose interpolation shape within the same four-point series. They are not independent validation or experimental uncertainty.

### Process and reduced-model checks

Maximum unit/overall relative mass residual: 2.231371e-13.
Maximum enthalpy residual normalized by total shaft power: 1.618426e-06.
Largest omitted-axial-momentum/friction upper screen: 1.317500% (limit 2%).

The checker replays 1620 source-capacity comparisons and 1458 discrete policy selections. Independent upstream replay recorded: `True`.

[Scientific figure](operating_map.png) · [Compact numerical reference](operating_map_reference.json) · [Reproduction commands](../operating_map/README.md)

The standalone map adds no new Panel browser controls. The existing low-Re interface retains its previous browser verification; no new browser pass is claimed for this CLI/report extension. Reacting and non-isothermal additions remain deferred.
