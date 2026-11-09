> Historical record of commit f3b925b, superseded by the stock-component case.
> Its custom cylinder-intersection helper has been removed under the project
> OpenLB requirements. The numbers below apply only to that earlier implementation.

# Verified 3D hydraulic benchmark — 8 November 2026

The native OpenLB capillary case, Python coordination/exports, and Julia
permeability translation are implemented and tested on this Mac. The suite
passed its numerical, geometry, scaling, time-step, conservation, and export
checks. The inferred permeability drives a Julia Darcy column calculation.
The existing DWSIM compatibility case has not yet been connected to this output.

Setup, exact launch/test commands, assumptions, acceptance thresholds, and
the versioned data contract are in [hydraulic/README.md](../hydraulic/README.md).
The accepted result snapshot is
`results/hydraulic/runs/20261108T164130Z-9f02cde8/summary.json`.
The compact observation is preserved in `docs/hydraulic_reference.json`;
full raw outputs and failed exploratory runs remain in ignored `results/`.

## Physical and numerical reference

A 50 µm radius circular pore runs along +x for 400 µm, surrounded by a declared
200 µm × 200 µm representative solid cell. The cell area is 4e-8 m² and porosity
is π/16. Viscosity is 0.001 Pa s and reference density is 1000 kg/m³. These are
prescribed constant benchmark properties. Inlet and outlet gauge pressures are
0.1 and 0 Pa. The OpenLB solver starts from zero velocity; it does not impose
the analytical velocity or evaluate an analytical flow formula.

The independent steady, fully developed Newtonian solution is
`u(r)=Δp(R²−r²)/(4μL)`, `Q=πR⁴Δp/(8μL)`. With the declared cell area,
`k=πR⁴/(8A_cell)`. Thus:

| Analytical quantity | Value |
| --- | --- |
| Total volume flow | 6.13592315154e-13 m³/s |
| Cell-based permeability | 6.13592315154e-11 m² |
| Pore-area-based permeability, different convention | 3.125e-10 m² |
| Maximum axial velocity | 0.00015625 m/s |

OpenLB 1.9r0 uses a serial double-precision D3Q19 BGK solve, curved Bouzidi
no-slip walls, and regularized local-pressure boundaries. Physical wall
intersections use the exact circle geometry. Python reads SI quantities and
Julia calculates `μ L_taps Q_mid/(A_cell Δp_taps)`. The pressure taps are at
x=100 and 300 µm; no analytical permeability replaces that measurement.

## Observed grid results

| Cells across diameter | Inlet / outlet nodes | Inferred k (m²) | k error | Velocity-profile relative L2 error |
| --- | --- | --- | --- | --- |
| 16 | 193 / 193 | 6.04852374096e-11 | 1.424% | 1.156% |
| 24 | 437 / 437 | 6.09076342218e-11 | 0.736% | 0.599% |
| 40 | 1245 / 1245 | 6.12435863056e-11 | 0.188% | 0.154% |

The finest measured flow is 6.12298269908e-13 m³/s (0.211% below the reference).
Its interior tap pressure drop is 0.0499887667301 Pa across 0.0002 m, and endpoint pressures
are recovered within the 1e-9 Pa check. Both physical endpoint locations and
expected cross-section node counts are checked before solving. Additional
geometry-only cases N=20, 32, 40 passed with 305, 793, 1245 nodes per boundary.
This includes a non-binary mesh spacing at N=40, which had exposed a missing
outlet in the separate stock 2D compatibility example.

On the N=40 run:

- Relative inlet/outlet **mass** mismatch: 7.0523e-8.
- Interior-plane mass-flux spread: 3.1139e-8 relative.
- All-plane mass-flux spread, including end planes: 2.4164e-4, or 0.02416%.
- Interior mean-pressure nonlinearity: 1.8204e-8 relative to the tap drop.
- Measured maximum Mach number: 6.7582e-5.

Flux integration sums lattice-center velocities and densities with weight Δx².
It applies no analytical area/flow correction. LBM density changes slightly
with pressure; consequently volume-flow mismatch alone is not the mass-balance
test. The exported diagnostics distinguish mass and volume flow, and include
the larger end-plane variation. These integrated tests do not claim exact
per-cell residuals or full finite-volume conservation.

## Scaling and numerical sensitivity

At N=24, the same solver binary produced:

| Check | Observation |
| --- | --- |
| Double Δp, 0.1→0.2 Pa | Q ratio 1.999947882; k shift 0.00261% |
| Reverse Δp, 0.1→−0.1 Pa | Q ratio −1.000052085; positive k shift 0.00521% |
| Half Δt, τ=0.8→0.65 at fixed grid/viscosity | k changes 0.1328%, below the 0.5% check |
| Tighten steady tolerance, 1e-8→1e-10 | k changes 2.43e-11 relative |

Small pressure-scaling differences are consistent with the weakly compressible
LBM formulation; this is not an exact incompressible sparse-network solve.
The time-step check also changes relaxation time and its interaction with the
boundaries; it does not establish a formal temporal convergence order.

## Boundary checks and rejected attempts

The generic OpenLB indicator ray search uses a tolerance of 1e-4 **of a lattice
link**, not an absolute SI length. An earlier absolute-tolerance diagnosis was
incorrect and was corrected after reading the underlying search implementation.
Its interpolated-boundary setter switches to half-way bounce-back when a returned
intersection lies outside the link. On the final N=40 geometry, diagnostic
evaluation of the generic search returned 12,681 out-of-link estimates among
111,744 wall queries.

The case uses exact quadratic circle intersections, a unit-consistent link norm,
and a guard for endpoint roundoff within 1e-10 of a link. The guard treats that
endpoint just inside q=1; it does not fit the wall to the flow reference.
There were 13,484 endpoint guards and zero failed exact intersections on N=40.
The side-wall condition also includes pressure-boundary fluid planes at the
circular rim. Including those planes reduced the N=16 all-plane mass-flux spread
from 0.662% to 0.277%, while barely changing interior permeability. No OpenLB
library source was patched.

Rejected coarse/earlier boundary variants remain in `results/hydraulic/`; they
are superseded by the final 16/24/40 study. No acceptance threshold was relaxed.
The selected three grids show decreasing errors and meet the finest-grid limits;
arbitrary-resolution monotonicity and an asymptotic convergence order are not
claimed. Further geometries require their own boundary and convergence checks.

All seven final cases were solved afresh with the same final executable.
Separate export verification checked their hashes, int64 index materialization,
SI coordinate agreement, table shapes, and Julia recomputation of all derived
quantities. The optional reuse route was exercised during development and
requires matching executable, inputs and artifact hashes; the final run used
no cached solver cases.

## Julia upscaling and regression checks

The finest inferred k supplies a separate 0.02 m column, area 1e-4 m² and drop
5 Pa, retaining the benchmark's pressure gradient. Julia returns
**1.53108965764e-9 m³/s**, or **0.0918654 mL/min**. This demonstrates the same
idealized straight-pore medium repeated at larger scale. It is not an empirical
material/device validation. The original network app's 1000 Pa column default
remains unchanged.

Direct Julia tests passed **208 assertions**: the previous 191 plus 17 measured
inference/reference tests. The existing Python/Julia/Panel suite passed
**11/11 tests**, using normal macOS graphics access. Fourteen native invalid or
unconverged cases returned clean errors and exported no completed raw result.
No new Python or Julia dependency was added; both saved dependency locks remain
unchanged. The benchmark figure was rendered with Matplotlib and visually
inspected. This CLI milestone did not add or visually verify a capillary UI.

Tested versions remain Python 3.12.9, Julia 1.11.6, JuliaCall/PythonCall 0.9.37,
OpenLB 1.9r0, GCC 15.2.0, tinyxml2 11.0.0 and macOS zlib 1.2.12 on macOS 26.6.2
ARM64. A new case build took 19.57 s with the existing OpenLB core; N=40 took
226.43 s to reach 16,800 steps. These are individual local timings, not a
performance comparison. In the restricted terminal sandbox Julia initially
could not update its normal user compilation cache; authorized host access
resolved that restriction. No installation or numerical blocker remains for
the documented benchmark.

## Next executable step

Connect the measured property record to a Julia Darcy pressure-drop calculation
at a declared flow, then feed that pressure requirement into the existing real
DWSIM pump automation. Match fluid properties, units, representative area and
operating regime explicitly. Introduce a second pore radius with fixed cell area,
rerun OpenLB, and verify that the resulting pressure/power change follows the
measured resistance. Those connections and the second-geometry comparison are
still pending; no DWSIM process result in this milestone used OpenLB permeability.
