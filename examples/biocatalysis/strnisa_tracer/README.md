# Strniša micro-packed-bed inert tracer benchmark

This is a separate, **nonreacting** reference for the next enzyme reactor.
It uses the measured blank measurement loop and the two dimensional glucose
pulses from [Strniša et al., 2018](https://doi.org/10.1016/j.cej.2018.05.028),
Fig. 3b. A (23.9 mm) is the qualification/calibration boundary; C (81.0 mm)
is a held-out length. See the [scientific report](REPORT.md),
[frozen physical protocol](protocol.toml) and
[source extraction](../../../literature/sources/strnisa2018/README.md).
No reaction rate or DWSIM process calculation is part of this test.

## Prerequisites

From the project root, use the established [native reproduction guide](../../../docs/guides/reproduction.md)
and Python/Julia environments. The executed native build uses OpenLB 1.9r0,
the stock D3Q19 TRT collision and forcing/bounce-back facilities, Homebrew
GCC 15 and local OpenMP. [The existing build helper](../../adsorption/transport_root_cause/hydraulics/stock_omp.py)
verifies vendor integrity and builds in `.local/biocatalysis/strnisa_tracer/`;
application sources remain in `openlb/`, outside the vendor tree.
The compact reference and source tables do not require the private article
PDF or any additional experimental data.

```sh
source scripts/env.sh
uv run --locked --offline pytest -q examples/biocatalysis/strnisa_tracer/test_loop.py
julia --threads=2 --startup-file=no --project=julia examples/biocatalysis/strnisa_tracer/test_transport.jl
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native build
STRNISA_RUN=$(uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native prepare)
```

The first H12 attempt used `tau_even=0.58` and failed the numerical density
variation screen. Retain its field and then warm-start the same geometry at
smaller lattice time step. It is important not to confuse this warm start
with a new accepted field. This failed attempt is expected on the tested
machine; if it completes on another machine, still keep it separately.

```sh
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native run "$STRNISA_RUN" A 12 --tau 0.58
# Inspect raw.json before preserving the rejected attempt (expected nonzero exit).
mv "$STRNISA_RUN/A_H12" "$STRNISA_RUN/A_H12_attempt_tau058"
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native run "$STRNISA_RUN" A 12 --tau 0.548 --warm "$STRNISA_RUN/A_H12_attempt_tau058"
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 12 handoff --diagnostic
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 12 trace --diagnostic --dt 0.00125
```

H12's local flow correction exceeded 5%, so it is retained as a diagnostic,
not a qualified tracer carrier. The independently passing H16 and H20
fields use a common physical marker time step. The short reactor's H16
half-step and diffusivity tests do not modify its saved flow.

```sh
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native run "$STRNISA_RUN" A 16 --tau 0.56
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 16 handoff
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 16 trace --dt 0.00125
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native run "$STRNISA_RUN" A 20 --tau 0.58
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 20 handoff
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 20 trace --dt 0.00125
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 16 trace --dt 0.000625
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 16 trace --dt 0.00125 --Dfactor 0.8
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" A 16 trace --dt 0.00125 --Dfactor 1.2
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.source_audit
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.compare "$STRNISA_RUN" freeze --depth 16 --dt 0.00125
```

The freeze checks H16–H20 mean/variance and the H16 half-step before
recording the complete short predictions, any A-only delay parameters, and
the source and code hashes. Its status distinguishes numerical qualification
from agreement with the complete experimental loop. If a screen fails on
another machine, investigate rather than overriding it silently.

```sh
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.native run "$STRNISA_RUN" C 16 --tau 0.54
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" C 16 handoff
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.handoff "$STRNISA_RUN" C 16 trace --dt 0.00125
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.compare "$STRNISA_RUN" predict
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.compare "$STRNISA_RUN" score
```

`predict` reads only the measured blank, the frozen A parameters and saved C
tracer travel times. `score` opens the measured C curve only after its
predictions have been saved. The independent tracer and related-bead storage
sensitivity are separate; none uses a reaction calibration.

## Saved contracts and checks

- `plan.json` specifies the two reactor dimensions, assumed geometry,
  porosity, physical fluid properties, seeds and frozen protocol hashes.
  Bead files contain `(x,y,z,radius)` in metres. The axial seam is open to
  the tracer, with developed flow obtained from the periodic pressure-driven
  full-length carrier; the physical lower plate is `z=0`.
- Native `field.f64` has little-endian Float64 rows `(i,j,k,u_x,u_y,u_z,rho)`
  with integer-valued zero-based lattice coordinates (lateral solid pad
  included), velocity in m/s and density in kg/m³. Flow is positive axial;
  pressure-gradient units are Pa/m. The saved native mask must agree exactly
  with the independent bead and plate rasterization. Native and marker time
  steps are different.
- Shared faces are oriented from the lower lattice coordinate toward its
  positive-axis neighbour; Julia node IDs are one-based in the face files.
  The rates are signed m³/s. The Julia handoff changes the rates only to
  preserve the native throughput and local balance; its correction norm is
  screened. Compact incoming/outgoing faces have shape `(3, N)` in Julia
  (row-major `(N, 3)` on the Python side), with zero neighbour ID for an
  impermeable wall or the open axial seam. The seam carries the saved
  forward flow and the no-diffusive-flux approximation.
- The inert marker files contain one little-endian Float64 exit age in
  seconds per injected marker. A non-arrival is recorded in the binary file
  as NaN, is **not** redistributed or counted as recovered, and is reported
  in the finite JSON mass ledger. The response is through the complete
  reactor, not a repeated cell. The blank and predictions use mM and s;
  the integral at 100 µL/min is the glucose amount in µmol. The normalized
  moments/CDF do not replace the absolute dose check.
- The measured blank's pulse already includes injection, tubing, detector
  and acquisition averaging. It is convolved once with the impulse
  response. The two extra-retention closures are A-calibrated diagnostics;
  their physical location is unqualified. The separate related-carrier
  internal-liquid model is a moment-only sensitivity, not a fitted pore law.

Full runs are ignored under `results/biocatalysis/strnisa_tracer/runs/`.
The exact run used by the report is identified in the small
[reference summary](reference/summary.json). To regenerate the tracked
compact evidence and figures from a completed new run, explicitly use:

```sh
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.reporting "$STRNISA_RUN"
uv run --locked --offline python -m examples.biocatalysis.strnisa_tracer.check --require-run
uv run --locked --offline python scripts/check_repository.py
```

Without the local full run, `check` (without `--require-run`) checks the
compact source, model and output hashes only; it does not claim to rerun
OpenLB or Julia.
