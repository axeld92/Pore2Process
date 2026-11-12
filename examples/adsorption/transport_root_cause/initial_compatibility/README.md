# Separate initial-condition compatibility check

The original sinusoidal concentration varies normally to a sphere at most
surface points. This is a valid weak initial/boundary-value problem and its
continuum solution must still respect the initial concentration range. It can,
however, create an immediate thin boundary adjustment that is difficult to
resolve. This control asks whether removing that particular initial adjustment
substantially reduces the numerical peaks.

Only initial concentration changes. Near each intended smooth sphere it retains
the original sine's surface trace and tangential variation, extends that value
constantly along the normal, and blends to the original field over a fixed
72.917333 micrometre buffer. The buffer is less than half the minimum particle
gap. Stock `IndicatorSphere3D`, stock vector normalization and stock
`PolynomialStartScale` provide the geometric and smoothing facilities. Geometry,
velocity, diffusivity, time step, collision and wall operators are unchanged.

The independent minimum-image FCC oracle agrees with all native initial N32
samples within 7.8e-16 kg/m³. Native surface probes verify wall trace,
periodicity, zero normal derivative and a smooth transition. All 11 construction
checks pass; five independent profile tests and one JSON typing test pass. Initial liquid concentration
is .250197933–.749802067 kg/m³ with mean .5 kg/m³.

Both tested grids still fail the concentration-range screen:

| Grid | Normal-compatible full-run C range, kg/m³ | Change in largest excursion from original sine |
| --- | --- | --- |
| N32 | .118326743–.881863418 | about 0.49% smaller |
| N64 | .125475664–.880145355 | about 1.25% larger |

The N64 buffer is the same 72.917333 micrometres, and all 11 initial-construction
checks pass there too. Normalized kinetic mass residuals remain below 7e-13.
Removing smooth-surface initial incompatibility alone does not materially
explain or remove the dominant error. The changed initial profile is a
diagnostic; its later concentration field need not match the original field.

The normal condition here concerns the **intended smooth sphere**. The actual
voxel/full-way BounceBack wall has different represented locations and
directions. Persistent errors therefore do not rule out initial-data effects
at that discrete wall, or continued wall-treatment and spatial-resolution
effects. A different initial condition cannot qualify the original sine case,
an open-boundary tracer or a column dispersion coefficient.

## Preserved execution and typed assessment

The first controller completed its native simulation successfully, then strict
JSON rejected NumPy boolean values in the assessment. Its original source,
binary, initial/final/peak fields and execution return code remain intact in
`runs/20261110T214547Z-e8d158cc`; `controller_failure.json` records that history.
The separate assessment converts NumPy scalar types to ordinary JSON values.
It changes no threshold, initial profile or simulation and does not rerun N32.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.initial_compatibility.assessment \
  --native-parent results/adsorption/transport_root_cause/initial_compatibility/runs/20261110T214547Z-e8d158cc \
  --grid 32
uv run python -m examples.adsorption.transport_root_cause.initial_compatibility.assessment \
  --check results/adsorption/transport_root_cause/initial_compatibility/assessments/20261110T214943Z-860fc2bb
uv run python -m examples.adsorption.transport_root_cause.initial_compatibility.assessment \
  --grid 64 \
  --native-parent results/adsorption/transport_root_cause/initial_compatibility/runs/20261110T214547Z-e8d158cc \
  --flow-summary results/adsorption/transport_root_cause/hydraulics/omp_runs/20261110T210317Z-fdf047ba/summary.json
uv run python -m examples.adsorption.transport_root_cause.initial_compatibility.assessment \
  --check results/adsorption/transport_root_cause/initial_compatibility/assessments/20261110T214946Z-5eab0ccf
uv run pytest tests/adsorption/test_initial_compatibility.py -q
```

The [protocol](protocol.md) declares the profile and comparison before execution.
`initial.csv`, `initial_profile_probes.csv` and metadata retain the actual native
initial data. `run.py` preserves the original controller; `assessment.py`
provides the typed continuation and explicitly guarded N64 execution using the
same archived binary and fixed physical buffer.

Exactly two native simulations completed. The N32 typed reassessment launches
none; N64 launches one new case. No native simulation was interrupted in this
component. The [compact reference](reference/compatibility_reference.json)
retains both source summaries, original comparison summaries and their hashes.
