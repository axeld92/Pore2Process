# 3D circular-capillary hydraulic benchmark

Python runs the native OpenLB 1.9 lattice-Boltzmann case, validates its SI exports
and sends measured flow and pressure drop to Julia. Julia infers axial Darcy
permeability. The separate [connected comparison](../coupling/README.md) uses
that property to drive actual DWSIM pump solves.

From `/Users/dula/Documents/ChatGPT/Pore2Process`, after the
[tool setup](../compatibility/README.md):

```sh
source scripts/env.sh
uv run python hydraulic/build.py
uv run python hydraulic/run.py
uv run python hydraulic/validate_inputs.py
uv run python hydraulic/check_exports.py
uv run python hydraulic/plot.py
julia --startup-file=no --project=julia julia/test/runtests.jl
uv run pytest -q
```

The check/plot commands default to the summary in `results/hydraulic/latest.json`.
They also accept an explicit `summary.json` path. Builds and results are ignored
in Git. Each suite keeps a unique run folder, including unsuccessful attempts.
`latest.json` reports the suite's status; inspect it and the command's exit code.
`--reuse /absolute/path/summary.json` requires matching binary, inputs and raw
artifact hashes, copies those artifacts and reruns Julia analysis. The default
makes fresh solves. No package installation happens on launch.

## Implementation boundary

All application code, headers and configuration are in `hydraulic/openlb/`.
The main follows the stock `poiseuille3d` example's mesh → geometry → lattice →
initial values → simulation → results lifecycle. Dynamics, velocity set, units,
geometry, boundary conditions, plane integration, and convergence monitors are
existing OpenLB components. There are no project numerical templates, custom
cylinder indicators or boundary algorithms. See [the policy and API map](../docs/openlb_code_policy.md).

The original OpenLB source is never written. Original files are checked against
the pinned official archive before/after the external build. Objects, executable,
core archive and configuration remain in `.local/hydraulic/build/`. This Mac
uses existing Homebrew GCC 15.2, tinyxml2 11 and system zlib. Paths in the external
configuration need deliberate adaptation on another machine.

## Physical model and independent reference

The default circular pore has radius 50 µm, length 400 µm, axis +x, centerline
at y=z=50 µm. The representative square cell has fixed width 200 µm and area
**4e-8 m²**; porosity is π/16. Uniform benchmark properties are μ=0.001 Pa s,
ρ=1000 kg/m³, with inlet/outlet gauge pressure 0.1/0 Pa. These prescribed values
are distinct from the actual DWSIM water properties used in the connected run.

The current solve uses serial CPU_SISD, double precision, stock D3Q19 TRT
(magic parameter 3/16), a stock cylinder
indicator and Bouzidi no-slip wall, including its rim at the pressure planes.
Flat local-pressure boundaries idealize fully developed flow. They do not model
reservoirs, physical entrance loss, gravity, chemistry or transient material
storage. Initial velocity is zero, and initial pressure is linear. No analytical
flow/velocity is imposed. LBM is weakly compressible; mass and volume diagnostics
are distinct.

The earlier BGK references are preserved as historical results. The TRT
regression is saved separately under `results/hydraulic_trt/`; all eight solver
cases, three geometry audits and the independent SI/export replay pass.
Reproduce that suite with `uv run python hydraulic/run.py --output results/hydraulic_trt`.
See [the TRT update](../docs/star_channel_showcase.md) for its numerical context.

For steady axial Newtonian flow, integrating
`μ(1/r)d/dr(r du/dr)=−Δp/L` with finite center velocity and u(R)=0 gives
`u(r)=Δp(R²−r²)/(4μL)`; integrating over the pore gives
`Q=πR⁴Δp/(8μL)`. Using Darcy superficial velocity Q/A_cell gives
`k_exact=πR⁴/(8A_cell)`. Using pore area instead gives R²/8, a different
permeability convention. These formulas are independent post-processing oracles.

Julia's `infer_permeability(L_taps,A_cell,μ,Δp_taps,Q_mid)` uses the actual
measured flow and pressure taps at L/4 and 3L/4. The default inferred property
also supplies a Julia column of length 0.02 m, area 1e-4 m² and pressure drop
5 Pa. This assumes the same idealized aligned medium repeats at larger scale.
The original Panel network example keeps its own 1000 Pa column default.

## Mesh and stop checks

The suite reports N=16,24,32,40 cells across the diameter. At fixed N and ν,
Δt=(τ−0.5)Δx²/(3ν); τ=0.8 is default and τ=0.65 halves Δt.
OpenLB's two `util::ValueTracer` monitors require relative sample standard
deviation below 1e-8 over six flow/pressure samples. Samples are separated by
max(100,N²/2) steps. A separate run uses tolerance 1e-10. Python independently
checks the recorded six-sample window. The 60000-step limit fails without
exporting a completed `raw.json`.

Circular-wall discretization and quadrature can produce cancellation on an
individual coarse grid. All four values are reported; monotonic error and a
formal convergence order are not claimed. Acceptance requires N32→40 k change
<0.1%, N24→40 change <0.5%, finest analytical k/profile errors <1%, and finest
all-plane mass spread <0.1%. The existing accuracy and conservation thresholds
are retained for the revised stock-component case.

All cases require finite results, requested input/pressure snapshots, typed
indices, SI coordinate agreement, valid shapes, relative inlet/outlet mass
mismatch ≤2e-5, interior-plane mass spread ≤1e-4, all-plane mass spread ≤0.5%,
and pressure nonlinearity ≤0.2%. Pressure doubling/reversal must change k by
<0.1%; halving Δt by <0.5%; tightening the stop tolerance by <1e-7 relative.
Extra geometry-only runs at N20,32,40 check counts and physical endpoint planes.
Twenty-five invalid/unconverged native cases must fail cleanly without a
completed result.

## Export contract, version 2

Pressure is gauge Pa; all dimensional fields are SI. Positive Q/mass flux is
+x on every plane, including inlet and outlet. The serial case uses one cuboid.
CSV indices are zero-based int64, physical fields float64, with
**x=iΔx, y=(j+0.5)Δx, z=(k+0.5)Δx**. Half-cell transverse centers avoid points
exactly on the circular wall. Rows run with j outermost, k innermost.
Version 1 used y=jΔx,z=kΔx and the earlier custom wall helper; its historical
results are preserved separately. Python can still read those original snapshots.

| Output | Contents |
| --- | --- |
| `raw.json` | Actual settings, dimensions, explicit axial interval count, transverse offset, wall/stop method, measured flow/pressure, time step, convergence and timing |
| `planes.csv` | Structured shape `(axial_intervals+1,)`: index, x, node count, Q, mass flow, mean pressure |
| `midplane.csv` | One row per fluid center at x=L/2: i/j/k, x/y/z, velocity, gauge pressure and density |
| `history.csv` | Steps/time, inlet/mid/outlet flow, inlet/outlet mass, tap pressure, successive changes |
| Per-case `summary.json` | Measurements, Julia-derived quantities, independent comparisons, command, reuse origin and raw hashes |
| Suite `summary.json` | Cases, acceptance, geometry, versions, source/binary/build hashes, column inputs |
| `analysis_verification.json` | Separate artifact, coordinate/index and Julia recomputation check |
| `hydraulic_benchmark.png` | Numerical/analytical velocity, pressure and four-grid error comparisons |

OpenLB's own plane-integral functors perform discrete integration over the
existing circle indicator. There is no custom quadrature or analytical area/flow
correction. Integrated conservation includes discretization, quadrature and
end-boundary effects; it is not an exact per-cell finite-volume residual.

Native application options include `--resolution`, `--radius-m`,
`--viscosity-Pa-s`, `--density-kg-m3`, `--pressure-drop-Pa`, `--tau`,
`--steady-tolerance`, `--max-steps`, and `--geometry-only`. Require even N8..64,
radius 20..60 µm, positive finite fluid properties, nonzero |Δp|≤1 Pa,
τ in (0.5,1.5], tolerance in (0,1e-3], and length/Δx an integer multiple of four.
Only the documented geometries/settings have passed the scientific acceptance;
input validation alone does not validate another flow regime or material.
