# Critical-flux-constrained crossflow benchmark

Julia reproduces the published García table, classifies a fixed target flux,
and adds a separate monotonic constraint from Chiu's figure-reading intervals.
Python coordinates actual DWSIM recycle cases using the existing stock-unit
adapter. Previously verified OpenLB support measurements supply membrane
pressure loads through the saved Julia translation. No new OpenLB case or
experimental fit is introduced.

## Run and check

From the project root on the tested Mac:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.filtration.literature_benchmark.run --require-sources
uv run python -m examples.filtration.literature_benchmark.check
uv run pytest -q
/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
```

The same Python 3.12.9, Julia 1.11.6, JuliaCall/PythonCall 0.9.37, DWSIM
10.2.15.0 and private .NET SDK 10.0.403 are reused. No dependency or lockfile
changes are needed. Julia's normal user compilation cache must be writable;
the VTK regression tests need normal host graphics access. The initial
restricted preflight's Julia import failure was a cache-permission diagnostic;
the host-access inventory and actual Julia runs distinguish it from a missing
runtime.

The default parent is `results/crossflow_recycle/latest.json`. On a clean
checkout, first reproduce the [graded-support](../graded_support/README.md),
[filtration](../once_through/README.md) and
[recycle](../crossflow_recycle/README.md) ancestors. All ancestors are verified
before new process solves. Optional `--parent path/to/summary.json` and
`--output results/another_benchmark` select different saved inputs/output roots.

Every execution creates a fresh run directory and retains failures. The
`latest.json` pointer changes only after the new study verifies. Omitting
`--require-sources` permits replay from the tracked article extraction when
the ignored inbox is absent; present sources must match their recorded hashes.
No article, package or runtime is downloaded during a benchmark run.

## Scope and data contract

[protocol.json](protocol.json) freezes the scenario choices. The
[source review](../../../literature/crossflow_papers_review.md) and
[extraction](../../../literature/crossflow_extraction.json) distinguish transcribed
tables, pixel readings and missing experimental conditions.

- Julia receives SI pressures, volumetric flows, areas, fluxes and concentrations.
  Python converts source LMH and L/h by dividing by 3,600,000. Row/point IDs
  remain one-based literature identifiers. Optional Julia source indices are
  one-based and are distinct from the original network's zero-based node IDs.
- The source energy is `delta_p*Q/(J_critical*S)`. Scenario circulation energy
  uses fixed target flux; the DWSIM result includes both pump shaft duties.
  Cooling heat is separate. No motor/refrigeration electricity is inferred.
- Flow factors **1 and 7** represent competing interpretations of the source
  Q column, not accepted corrections. Axial loss stays fixed. Parallel area
  normalization scales flow by process/source area without altering axial loss.
  The resulting pointwise `delta_p/Q` values are not a fitted linear channel law.
- Decisions use measured table operating points. The separate figure policy
  assumes nondecreasing critical flux and takes only definitely-left/right
  measurements as bounds, accounting for the recorded reading intervals.
  It does not interpolate through the low-speed transition or extrapolate.
  `indeterminate` means insufficient evidence for that policy, not infeasible.
- Capacity fractions 1.0/0.8 are explicit policy sensitivities. Lowest energy
  means lowest among solved admissible discrete points, not a continuous optimum.
- The unchanged 100 LMH, 1 m², 50% recovery, 298.15 K water-carrier scenario
  reuses native support resistance and the fixed literature-informed non-support
  impedance. The same empirical filtration criterion applies to each support
  design. Geometry-dependent fouling has not been validated.

Julia also balances a passive solute through the finite-pass membrane,
retentate split, mixer and whole loop. Module-inlet concentration is held at
0.01 kg/m³ with assumed rejection 1. Fresh concentration is adjusted to match
that reference condition as circulation changes; this is not a fixed-fresh-feed
composition sweep. The model explicitly exports all concentrations and three
component residuals. DWSIM still contains Water only. Suspension rheology,
particle transport and a real solute property/rejection model remain outside
this conditional dilute-carrier demonstration.

## Numerical and simulator checks

The matrix has 30 design cases (five nominally admissible source rows, two
flow conventions and three preserved support designs) plus two high-circulation
seed/tolerance checks. A separately forced one-iteration failure must be rejected.
All DWSIM cases save, reload and solve using stock pumps, mixer, cooler, recycle,
splitters and pressure-loss valves. Balance, property, convergence and archive
connection checks reuse the established validator.

The factor-seven scenario reaches a module/permeate ratio of approximately
187. Stock direct substitution then needs more iterations than the previous
prototype's limit. Only our project adapter's input ceiling was raised from
1000 to 5000; its existing default and stock algorithm are unchanged. The new
protocol requests 5000 and tests a low-flow seed and a high-flow seed with
tenfold tighter tolerances at this largest ratio. Original OpenLB and DWSIM
distribution files are untouched.

Exports include immutable input/extraction/ancestor snapshots, source and
implementation hashes, actual saved DWSIM flowsheets and renders, stock-unit
states, source admissibility decisions, passive-solute balances, comparisons,
discrete selections, runtime versions and cold-initialization timing. Strict
JSON uses null for unavailable one-sided flux bounds. The checker replays
Julia and inspects all saved evidence; it does not launch external solvers.

See [the result and limitations](../../../docs/crossflow_benchmark.md). The original
Panel capillary-network example remains separately launchable.
