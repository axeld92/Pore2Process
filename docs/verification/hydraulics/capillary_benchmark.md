# Stock-component OpenLB hydraulic benchmark — 8 November 2026

> **Dated verification record.** This retains the execution/audit history.
> Use the [maintained scientific report](../../../examples/hydraulics/capillary_benchmark/REPORT.md) for the current explanation
> and the example README for current commands.

## TRT revision — 9 November 2026

Porous-media runs now require at least TRT. The current source uses stock
`TRTdynamics` with magic parameter 3/16. Run
`20261109T145749Z-e4c8d041` under `results/hydraulic_trt/` passes all eight
solver cases and three geometry audits. Its separate export checker verifies
recorded artifact hashes, SI coordinates/indexing and real Julia recomputation.
The original numerical screens are unchanged.

The finest N40 permeability is 6.13304992846e-11 m² (0.04683% analytical error),
with 0.04990% velocity-profile L2 error and 0.03361% all-plane mass spread.
The N32-to-N40 permeability change is 0.02052%. Pressure reversal/doubling,
the smaller time step and tighter steady criterion all pass. Cold Julia
initialization was recorded separately (2.77 s). Source/binary hashes identify
the TRT implementation even though this inherited raw-export schema has no
separate collision-name field. These values are solver verification results,
not a speed benchmark or material validation.

The sections below preserve the earlier BGK run and its historical numbers.
The new [star-channel showcase](../filtration/star_channel.md) also revalidates the
porous supports with explicit TRT metadata before new process comparisons.

## Earlier stock-component result

The revised 3D capillary case passes its analytical, mesh, conservation,
pressure-scaling, time-step, steady-state, geometry and export checks. It uses
OpenLB's existing numerical components throughout and follows the stock example
lifecycle. The custom cylinder helper is removed. All additions/build products
stay outside the original distribution; 2,070 original files match the official
archive. See [the code policy](../../architecture/openlb_policy.md) and
[launch commands/data contract](../../../examples/hydraulics/capillary_benchmark/README.md).

The previous implementation and measurements remain in Git and
[the historical report](capillary_legacy.md). Its numbers apply to that
implementation only. Earlier accepted run: `20261108T181448Z-84563000`.
Full exports: `results/hydraulic/runs/20261108T181448Z-84563000/summary.json`; compact observation:
`examples/hydraulics/capillary_benchmark/reference/hydraulic_reference.json`.

## Model and observations

Radius 50 µm, pore length 400 µm, fixed cell area 4e-8 m², viscosity 0.001 Pa s,
reference density 1000 kg/m³, and gauge pressure 0.1→0 Pa. The independent
Hagen–Poiseuille oracle gives Q=πR⁴Δp/(8μL)=6.13592315154e-13 m³/s and
k=πR⁴/(8A_cell)=6.13592315154e-11 m². The actual numerical solve starts at rest;
it does not impose the analytical velocity or use analytical k in translation.
Julia infers k from measured midplane flux and interior pressure taps.

| N across diameter | Nodes per boundary | Inferred k (m²) | k error | Profile relative L2 error |
| --- | --- | --- | --- | --- |
| 16 | 208 | 6.13161000918e-11 | 0.0703% | 0.1837% |
| 24 | 448 | 6.10808527003e-11 | 0.4537% | 0.4168% |
| 32 | 812 | 6.13090118395e-11 | 0.0818% | 0.0971% |
| 40 | 1264 | 6.13088508149e-11 | 0.0821% | 0.0781% |

N32→40 changes k by 0.000263%;
N24→40 by 0.3733%. Errors are non-monotonic:
the N16 net permeability error benefits from cancellation, and N32/N40 are very
close. We report every grid, check bounded fine-grid variation and the
independent oracle, and claim no formal convergence order. The original 1%
finest accuracy and conservation limits remain; monotonic net error is not a
valid assumption for this revised discretization.

Finest flow: 6.13078401407e-13 m³/s; measured tap drop:
0.0499991757518 Pa across 0.0002 m.
The stock plane-integral functors integrate physical field values without an
analytical area correction. N40 diagnostics:

- Inlet/outlet relative mass mismatch: 8.87333e-08.
- Interior-plane mass spread: 4.20784e-08 relative.
- All-plane mass spread, including endpoints: 0.03258%.
- Interior pressure nonlinearity: 5.32232e-09 relative.
- Maximum Mach number: 6.75424e-05; maximum-velocity Re: 0.0155982.

These are integrated checks, not exact per-cell residuals. Weak compressibility
means density-weighted mass and volume fluxes need separate diagnostics.

## Sensitivity and implementation checks

N24 baseline k=6.10808527003e-11 m²:

| Changed setting | Inferred k (m²) | Relative k change from N24 baseline |
| --- | --- | --- |
| double_pressure | 6.10792609846e-11 | 0.002606% |
| reverse_pressure | 6.10840340775e-11 | 0.005208% |
| half_time_step | 6.10068807239e-11 | 0.121105% |
| tighter_steady | 6.10808527018e-11 | 0.000000% |

Pressure doubling/reversal meet the 0.1% k-change check; the half-time-step case
meets 0.5%; the tighter stop case meets 1e-7 relative. Two stock `ValueTracer`
monitors check flow and pressure relative standard deviation over six samples;
exported history independently reproduces that stopping criterion.

The half-cell transverse grid x=iΔx, y=(j+0.5)Δx, z=(k+0.5)Δx avoids exact
node/wall coincidences while using the unmodified cylinder distance routine.
No replacement dynamics, descriptors, numerical templates, indicators, boundary
processors, plane integrators or convergence classes are present.

Eight fresh solver cases and three geometry-only cases passed. Separate export
verification checks hashes, SI coordinates, array shapes, int64 indexing and
Julia recomputation. Twenty-five invalid/unconverged native cases fail cleanly.
Direct Julia tests pass **215 assertions**, and the original Python suite passes
**11/11** with normal macOS graphics access. The original app, bridge, dependency
locks and Julia manifests are unchanged. The numerical figure was rendered and
visually inspected. This CLI milestone adds no capillary UI.

Python 3.12.9, Julia 1.11.6, JuliaCall/PythonCall 0.9.37, OpenLB 1.9r0, GCC15.2,
tinyxml2 11, zlib1.2.12, macOS26.6.2 ARM64. N40 used 16800 steps;
observed solve/export time 232.57 s. Runs overlapped on this Mac;
these timings are not a performance comparison. Julia initialization is recorded
separately in each suite summary.

The finest k supplies a 0.02 m Julia column, area 1e-4 m² and Δp=5 Pa:
Q=1.53272127037e-09 m³/s. The separate
[connected comparison](darcy_pump.md) now uses two OpenLB geometries
and actual DWSIM water properties to drive real pump calculations.
