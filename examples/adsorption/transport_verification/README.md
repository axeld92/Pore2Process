# Separate bulk-stretching and curved-wall checks

This continuation of the [packed-bed diagnosis](../transport_root_cause/REPORT.md)
separates three-dimensional deformation of a concentration field from transport
next to a single curved impermeable surface. The [prospective protocol](protocol.md)
defines the physical problem, the independent references and the screens.
The velocities are prescribed. No new bed, hydraulic solve, adsorption fit or
DWSIM simulation is implied by either test.

Native code follows the installed `advectionDiffusionReaction/advectionDiffusion1d`
and `laminar/poiseuille3d` examples. It uses existing OpenLB 1.9r0 geometry,
descriptors, scalar-specific regularization with stock TRT, and BounceBack.
Only the physical cases and their diagnostics are added here. The build stays
under `.local/adsorption/transport_verification`; the original distribution is
verified against the pinned archive. Python's continuum calculations are
verification references. They neither update OpenLB nor supply a closure.

## Run and replay

Tested environment: macOS ARM64, Python 3.12.9, NumPy 2.5.4, Homebrew GCC 15.2,
OpenLB 1.9r0 CPU_SISD/OpenMP. Uses the existing locked project environment.
No SciPy, environment changes or external data are needed.

```sh
source scripts/env.sh
OPENBLAS_NUM_THREADS=2 uv run python -m examples.adsorption.transport_verification.run --stage both
uv run pytest tests/adsorption/test_transport_verification.py -q
```

The first stage checks the independent reference, builds the stock executable
externally, then runs five native cases. The second stage uses the exact archived
executable from the first and runs eight cases. Each stage writes its own
timestamped folder under `results/adsorption/transport_verification/`. A failed
scientific accuracy screen does not erase a completed result. A failed input,
reference, structural or mass check stops the stage. If the first stage is
already complete, run only the second with `--stage curved --stretching-summary
<stretching-summary.json>`. Use `--check <summary.json>` to verify sources, saved
fields and analysis without another OpenLB run.

## Saved data contract

Every case contains exact launch arguments, a native exit status, model
parameters, executable fingerprint, every-step extrema, fluid inventory and
active-wall inventory. `material.u8` has material IDs 1 (fluid) and 2 (solid).
`initial.f64`, `middle.f64`, `peak.f64` and `final.f64` store concentration in
kg/m³. Arrays are headerless little-endian IEEE Float64 (except the unsigned
byte material mask), shape `(N,N,Nz)`, with z contiguous and zero-based cell
IDs. Velocity has shape `(N,N,Nz,3)`, x/y/z components in m/s. Cell centres
are `(index+0.5)*dx`; `raw.json` stores dimensions, spacing and time. The middle
snapshot is at `floor(steps/2)*dt`; the peak is the worst transient excursion
outside the prescribed range, or the initial state if there was no excursion.
The C values stored in solids are numerical populations and are excluded from
all physical concentration errors.

`stencil.csv` exports the installed OpenLB directions and weights. These are
used to verify the curved case's active reflected links and initial numerical
storage independently of the native loop. `reference_modes.npy` retains the
complex stretching reference at the two observation times. The disk oracle
is analytical, evaluated at the native fluid centres without fitting an
effective radius. The 50 µm near-wall shell is fixed in physical units.

Runtime evidence is intentionally kept out of Git; the report and compact
reference summary retain the scientific results and the exact run identities.
See the completed [scientific report](REPORT.md) and
[compact result record](reference_summary.json). Exact replay and plot commands
for the saved runs are:

```sh
OPENBLAS_NUM_THREADS=2 uv run python -m examples.adsorption.transport_verification.run \
  --check results/adsorption/transport_verification/stretching/20261111T004319Z-61b04382/summary.json
OPENBLAS_NUM_THREADS=2 uv run python -m examples.adsorption.transport_verification.run \
  --check results/adsorption/transport_verification/curved/20261111T011153Z-1a03a0b9/summary.json
OPENBLAS_NUM_THREADS=2 uv run python -m examples.adsorption.transport_verification.plot
```
