# Verified OpenLB → Python/Julia → DWSIM demonstration

On 8 November 2026 this Mac completed a reproducible end-to-end hydraulic
comparison with all four required tools. Changing pore radius changes the
measured permeability, Julia pressure requirement and actual DWSIM pump power.
This meets the first connected technical demonstration criterion for an
idealized synthetic medium. The subsequent application is now
[fixed-flux ceramic filtration](filtration_research_case.md); empirical material
validation remains unfinished.

Run ID: `20261108T181450Z-e1b190f9`. Full snapshot:
`results/coupling/runs/20261108T181450Z-e1b190f9/summary.json`. The compact observation is
`docs/coupling_reference.json`. [Exact commands, contracts and limitations](../examples/hydraulics/darcy_pump/README.md).

## Observed process effect

Both designs use a fixed 200 µm square representative cell and a 400 µm pore
length. The upscaled slab is 1 m long, area 1 m², at a fixed water flow of
1e-5 m³/s (0.6 L/min) and pump efficiency 75%.

| Pore radius (µm) | Measured permeability (m²) | Julia porous pressure load (Pa) | Actual DWSIM shaft power (mW) |
| --- | --- | --- | --- |
| 50 | 6.130870215e-11 | 145.18045 | 1.9357394 |
| 40 | 2.511215119e-11 | 354.44296 | 4.7259061 |

The smaller-pore permeability ratio is **0.409601742**;
the actual DWSIM power ratio is **2.44139587**.
Independent R⁴ theory predicts 0.4096 and 2.44140625. The measured property,
rather than that formula, supplies every downstream calculation.

Python launches actual OpenLB solves for both geometries at N24 and N40.
It passes measured signed flow and tap pressures to the independent Julia
package. Julia computes k=μL_taps Q/(A_cell Δp_taps), then Δp=μL_column Q/(kA_column).
Python sets the actual DWSIM pump's absolute discharge pressure to 100000 Pa+Δp,
and mass flow to ρQ. DWSIM solves its pump, material and energy streams, saves
the `.dwxmz` flowsheet, reloads and solves it again. Pressure and mass flow agree
after reload. The porous load is calculated externally; the flowsheet contains
a pump, not a new porous-medium unit operation.

DWSIM's Steam Tables (IAPWS-IF97) at 298.15 K and 100000 Pa supplies
μ=0.0008900825175722 Pa s and
ρ=997.0474354081 kg/m³. Those same values configure
OpenLB, Julia and mass-flow conversion. Gauge and absolute pressure conventions,
SI units and efficiency conversion are explicit.

## Validation and practical limits

- All four OpenLB geometry/grid solves pass analytical profile/permeability,
  integrated conservation, coordinate/indexing, finite-value and stopping checks.
  Both fine-grid errors are below 1%, with smaller k error than N24.
- The measured permeability ratio agrees with R⁴ within 1%; actual DWSIM power
  agrees with QΔp/η and its ratio with inverse measured k within 1e-9 relative.
- Mass balance is checked; liquid state, matching inlet properties and the pump
  energy stream pass. Actual energy residuals are
  7.92349e-07 and 7.92352e-07 W,
  or 0.0409% and
  0.0168% of shaft power.
  The 0.1% criterion accounts for the small enthalpy increment in the water flash;
  exact energy closure is not claimed.
- Estimated process maximum pore Re is 0.01141 and
  0.014263, below 0.05. Actual outlet viscosity/density
  changes remain below 0.01%, supporting the constant-inlet-property translation.
- Raw exports, actual inputs, translation snapshots, settings, diagnostics,
  timings, versions, source/binary/runtime hashes and saved flowsheets are kept.
  Separate verification recomputes the Julia translation and checks artifact hashes;
  a deliberately modified copy of an export is rejected.
- The OpenLB implementation follows the official example structure and reuses
  stock numerical components. Project code/configuration/build output remain
  separate; original source files match the pinned archive. See [the policy](openlb_code_policy.md).

The scientific comparison figure was visually inspected. This milestone is a
CLI workflow; no new connected UI or DWSIM desktop interaction is claimed.
The existing Panel network app remains intact, including its previously
recorded outstanding manual rotation-gesture check. No installation/coupling
blocker remains on this Mac. The medium is idealized straight, aligned pores;
there is no empirical material validation, isotropic permeability claim,
reservoir/entrance loss model, or full equipment model for a porous device.

## Next executable step

The selected question is now support geometry versus filtration pumping
energy at fixed flux and selective layer. A synthetic layers-in-series study
has passed; see [its results and next step](filtration_research_case.md).
The next geometric extension should vary
one geometric feature such as a constriction or tortuosity, validate the
OpenLB resolution/domain/pressure range, and propagate its measured resistance
through this same saved-data boundary. Decide whether that application needs
a pressure-loss unit in DWSIM. Cantera is optional for a later chemistry case;
it is not required for this hydraulic demonstration.
