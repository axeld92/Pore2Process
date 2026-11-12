# Native packed-bed tracer to Julia contact benchmark

[Scientific report](REPORT.md) · [Declared protocol and failed trials](protocol.md).

This case transfers actual OpenLB inert-tracer exit-flux information into an
independent Julia contact model. It does not change the batch-qualified particle
model, calibrate either Alhamed column, or launch a new DWSIM case.

## Prerequisites and commands

Use the existing macOS ARM64 setup and shared locked Python/Julia environments.
The OpenLB 1.9r0 archive/distribution and Homebrew GCC15 are needed by the shared
external build helper. The accepted FCC flow fields at24/32 and their bed summary
must exist. On a fresh machine follow the [parent run guide](../README.md),
including `uv run python -m examples.adsorption.bed`, first. Existing accepted
runs and vendor files must not be overwritten to force a different request.

From the project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/preflight.py --output results/adsorption/packed_bed_tracer/preflight.json
julia --project=julia julia/adsorption/test_tracer_contact.jl
uv run python -m examples.adsorption.packed_bed_tracer.run
uv run python -m examples.adsorption.packed_bed_tracer.check
uv run python -m examples.adsorption.packed_bed_tracer.conservative
uv run python -m examples.adsorption.packed_bed_tracer.conservative --check
uv run python -m examples.adsorption.packed_bed_tracer.plot
uv run python -m examples.adsorption.packed_bed_tracer.plot_conservative
uv run pytest tests/adsorption/test_packed_bed_tracer.py tests/adsorption/test_conservative_tracer.py -q
uv run python scripts/check_repository.py
```

`run` creates a unique folder, builds separate unmodified OpenLB core/application
objects in `.local/adsorption/packed_bed_tracer*`, executes one new N=16 flow and
eight scalar cases, then calls Julia on the main thread. Two independent native
subprocesses may run concurrently. JuliaCall never runs on those worker threads.

`check` verifies hashes and saved native balances/fields, and replays Julia
translation. It launches zero OpenLB simulations and does not replace saved
native evidence. Full native re-execution requires another `run` invocation.
`plot` makes figures from saved outputs. The existing full UI regression suite
requires normal macOS host graphics access; this case adds no browser interface.

## Scientific and file boundary

Flow is along positive x; y/z are periodic. Scalar reservoir planes enclose an
interior control volume; its link-surface separation is `(nx-2)*dx`, while the
reservoir-node separation is `(nx-1)*dx`. Physical positions in `final.csv` use
OpenLB coordinates: the first reservoir node is at x=dx/2. Integer i,j,k are
zero-based. Radius and packing are synthetic FCC, nominal porosity0.6. The finite
scalar segment reuses periodic fully developed pore velocities, including its
end-plane pattern; hydrodynamic entrances and collectors are not resolved.

- `request.json`: exact command, native input/source/binary hashes and SI inputs.
- `raw.json`: descriptor, collision, relaxation channels, dx/dt, physical domain,
  actual flow, interior fluid volume and timing. `length_m` is reservoir-node
  separation; translated `transport_length_m` is the control-surface separation.
- `flux.csv`: midpoint time(s), net inlet/outlet scalar flux(kg/s), positive in+x.
  These include diffusion and need not equal Q*C at the absorbing outlet.
- `inventory.csv`: fluid tracer, temporary wall-population tracer, cumulative
  feed/exit mass, raw conservation residual and concentration extrema.
- `final.csv`: i,j,k,material,coordinates(m), kinetic concentration(kg/m3),
  native physical pore velocities(m/s). Wall populations are numerical storage,
  not adsorption. Boundary populations differ from prescribed Dirichlet momenta.
- `julia_*.json`: explicit ADE SI inputs, trajectories and comparison metrics.
- `summary.json`: case diagnostics, numerical changes, contact derivation,
  versions, source snapshots/hashes and distinction between verification stages.

The normalized signal is exit flux divided by its late steady plateau. Its
moments are a boundary-specific transport response. Diffusion-dominated exit
transit is not interchangeable with fluid void residence time. Source Langmuir
parameters, particle diffusion, active inventory and experimental bed density
are never fitted or substituted here. Missing or failed qualifications remain
explicit; higher-Peclet exploratory runs are not silently promoted.

Direct Julia use, without Python or the interface:

```julia
include("julia/adsorption/TracerContact.jl")
using .TracerContact
absorbing_moments(length_m=0.004, velocity_m_s=0.0005, dispersion_m2_s=1e-6)
solve_ade(length_m=0.004, velocity_m_s=0.0005, dispersion_m2_s=1e-6,
          t_end_s=35., dt_s=0.002, cells=800)
```

The selected `conservative` stage holds measured Q and physical pore volume
fixed and derives only D from the short tracer mean. Its independent --check
replays six Julia predictions; no native or DWSIM solve is repeated. The older
two-moment v/D route remains a diagnostic shape comparison. The selected route
also tests absolute exit flux, not only separately normalized curves.
