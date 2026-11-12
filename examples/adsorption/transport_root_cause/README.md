# Transport resolution and root-cause investigation

This is a separate diagnostic continuation of
[transport accuracy](../transport_accuracy/README.md). It preserves the completed
examples and their failed data. The [prospective protocol](protocol.md) defines
the new spatial series, peak snapshots and exact local population accounting.
The [scientific report](REPORT.md) explains the engineering meaning and keeps
the failed particle-bed qualification explicit.

The particle bed is synthetic FCC geometry, not a reconstruction of the paper's
carbon. The physical application-scale flow is unchanged. D=1e-9 m2/s is the
declared molecular-diffusion scenario; no new experimental measurements or
experimental validation are claimed. The Julia particle model and DWSIM pumping
calculations are not changed or recalibrated by these tests.

## Components

- `study.py`, `openlb/`: matched-time-step spatial study, full peak fields and
  deterministic prefix replays. The native transport uses stock OpenLB
  `CombinedAdvectionDiffusionRLBdynamics` with scalar TRT and D3Q27, and stock
  full-way BounceBack.
- `parallel.py`: same scalar application using stock OpenMP, after a complete
  N32 serial/parallel field-equivalence check. Four threads are set in this
  process only. No numerical method is replaced.
- [hydraulics](hydraulics/README.md): independently converged fixed-flow fields
  and separately checked stock OpenMP execution.
- [periodic shear control](channel/README.md): isolates velocity variation
  without particles, compared with an independently converged continuum oracle.
- [fine shear controls](channel/thin/README.md): independently checked depth
  reduction, N96/N128 spatial refinement and a fine-grid time-step comparison.
- [flat-wall channel](flux_audit/flat_channel/README.md): separates smooth shear
  from impermeable wall treatment, with a qualified no-flux continuum reference.
- [population audit](flux_audit/README.md): actual local collision/streaming
  link balances; [localized findings](flux_audit/findings.md).
- [second-order candidate](flux_audit/candidate/README.md): tested stock
  equilibrium composition, rejected as a correction for the observed failures.
- [half-way reflection candidate](flux_audit/half_way/README.md): independently
  checks immediate stock reflection at the same voxel wall; rejected after
  larger porous excursions on N32 and N64.
- [N96 extension protocol](resolution96_protocol.md): conditional finer-grid
  particle-bed test after the wall-free N96 reference passes.
- [N96 time-check rule](resolution96_time_decision.md): recorded before the
  fine scalar execution; defines when a half-time-step confirmation is needed.
- [grid alignment](grid_phase/README.md): translates the same physical geometry
  and concentration to measure sensitivity to voxel representation.
- [initial compatibility](initial_compatibility/README.md): removes the initial
  smooth-sphere normal gradient without changing the transport operators.
- [matched background diagnostic](constant_protocol.md): uses the exact N64
  executable and time window to bound uniform leakage's contribution to the
  sinusoidal error; the diagnostic subtraction is never a model correction.

## Runtime and data contract

Target checked again: macOS 26.6.2 ARM64, Python 3.12.9, Julia 1.11.6,
NumPy 2.5.4, JuliaCall/PythonCall 0.9.37, GCC 15.2, unmodified OpenLB 1.9r0.
Existing `.venv`, `uv.lock`, Julia project/manifest and selected runtime are
reused. No dependency installation or lockfile modification is needed.
The actual inventory is saved at
`results/adsorption/transport_root_cause/preflight.json`.

Native fields are CSV with zero-based `i,j,k`, x-major/y/z ordering and SI
coordinates at `(index+0.5)*dx`. Material 2 is numerical solid; other material
labels are physical fluid in this all-periodic case. Concentration is kg/m3,
velocity m/s, length m, time s, inventory kg and signed scalar link flux kg/s.
Actual inputs and amplitudes are authoritative in `request.json`/`specifications`
and `peak.json`; the inherited textual `raw.json.initial_profile` string always
spells `.25` and is not the amplitude input for constant controls.

`inventory.csv` samples extrema and inventories every step. `peak.csv` is the
entire field at the time of the largest bound excursion, not the final field.
Whole-history minimum and maximum need not occur at precisely the same step.
`moments.csv` reports liquid-only mean and standard deviation. Its harmonic
coefficients are `2*mean_fluid[(C-0.5)*sin(2*pi*x/L)]` and the analogous cosine
projection. Because the fluid mask has a varying cross-section, their magnitude
is a weighted diagnostic, not a normalized Fourier amplitude or a fitted axial
dispersion. Numerical wall population inventory is separate from physical
liquid inventory.
Legacy `first_mode_amplitude` keys in the frozen grid-phase cohort refer to
that projection norm; they do not denote an ordinary Fourier amplitude.

The imposed superficial velocity corresponds to 43.8 mL/min through the 11 mm
reference column. The simulated periodic cell has transverse area `side_m²`:
its `raw.json.flow_m3_s` is the flow through that smaller area. In refined
source-point records, `actual_flow_m3_s` is the equivalent reference-column
flow. They must not be compared as if they had the same area. Native hydraulic
CSVs use the computational carrier (nu=1e-6 m2/s, rho=1000 kg/m3); the recorded
`velocity_scale` converts them to the application carrier at unchanged Reynolds
number before scalar initialization. Water population audits use the native
carrier's density/time/length scales; scalar velocities and times are physical
application SI values.

`audit_<step>.csv` records actual post-collision population links, predicted
streamed density and observed density. Shifted OpenLB populations are restored
before positivity checks. These diagnostics do not change the transport update.
`uniform_stencil.csv` predicts a first-order-equilibrium scalar's first step;
it is distinct from the actual water population flux audit.

## Execution

From the project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/preflight.py
uv run python -m examples.adsorption.transport_root_cause.study --grids 32 48
uv run python -m examples.adsorption.transport_root_cause.study --grids 32 48 --constant --end-s 3.803491159281669e-5
```

The saved source-point file for a finer mesh is only emitted after its actual
native hydraulic solution passes all fixed-flow, stationarity, mass, Mach,
density and porosity checks. Use it with its matching hydraulic summary:

```sh
uv run python -m examples.adsorption.transport_root_cause.parallel --grid 64 --point-json <qualified-run>/source_point.json --flow-summary <qualified-run>/summary.json
```

`<qualified-run>` is an explicit path under the ignored hydraulic results, not
a supplied dataset. See the hydraulic run guide for exact native setup commands.
The completed N96 execution used the following accepted continuation:

```sh
uv run python -m examples.adsorption.transport_root_cause.parallel --grid 96 \
  --point-json results/adsorption/transport_root_cause/hydraulics/continuations/20261110T223616Z-3f4c960a/source_point.json \
  --flow-summary results/adsorption/transport_root_cause/hydraulics/continuations/20261110T223616Z-3f4c960a/summary.json
uv run python -m examples.adsorption.transport_root_cause.parallel \
  --check results/adsorption/transport_root_cause/parallel/runs/20261110T230851Z-6bb12c48/summary.json
```

The complete serial/parallel scalar cross-check is:

```sh
uv run python -m examples.adsorption.transport_root_cause.parallel --grid 32 --verify-serial results/adsorption/transport_root_cause/runs/20261110T205424Z-ed34b563/summary.json
```

Prefix replays use the same executable and physical inputs, with
`--audit-step <peak-step> --end-s <peak-time-s>`; they are explicitly shorter
diagnoses and cannot qualify the full-duration case.

For a completed single-case parallel parent, `replay.py` also archives and
reuses its exact executable without rebuilding. Both modes retain the parent
hash, input snapshot and unchanged physical diffusivity:

```sh
uv run python -m examples.adsorption.transport_root_cause.replay --parent <parallel-run>/summary.json --peak-audit
uv run python -m examples.adsorption.transport_root_cause.replay --parent <parallel-run>/summary.json --half-dt
uv run python -m examples.adsorption.transport_root_cause.constant --parent <parallel-run>/summary.json
uv run python -m examples.adsorption.transport_root_cause.constant --check <constant-run>/summary.json
```

The peak replay explicitly checks that its final field reproduces the parent's
peak field. The half-step trial checks the same physical end time and reports
the direct physical-liquid field difference. It does not change the spatial
mesh, supplied velocity or molecular diffusion.

Saved-evidence checks and comparison (each summary argument is explicit):

```sh
uv run python -m examples.adsorption.transport_root_cause.study --check <scalar-run>/summary.json
uv run python -m examples.adsorption.transport_root_cause.parallel --check <parallel-run>/summary.json
uv run python -m examples.adsorption.transport_root_cause.assess <scalar-run>/summary.json <other-run>/summary.json
uv run python -m examples.adsorption.transport_root_cause.plot
uv run python -m examples.adsorption.transport_root_cause.check
uv run python -m examples.adsorption.transport_root_cause.assess --manifest
uv run pytest tests/adsorption/test_transport_peak.py tests/adsorption/test_transport_flux.py tests/adsorption/test_transport_channel.py tests/adsorption/test_transport_refinement.py tests/adsorption/test_transport_flat_channel.py tests/adsorption/test_transport_grid_phase.py tests/adsorption/test_initial_compatibility.py tests/adsorption/test_transport_accuracy.py tests/adsorption/test_transport_diagnosis.py -q
uv run pytest tests/adsorption/test_transport_half_way_review.py -q
```

These checks preserve and replay the archived evidence; they launch zero new
native simulations. The `study`/`parallel` commands without `--check` launch
actual simulations. Original OpenLB files, older numerical sources, dependencies
and unrelated files remain untouched. No UI/browser verification is implied by
this numerical diagnosis.
