# Column transport and hydraulic process qualification

Read [the scientific report](REPORT.md) and [the pre-evaluation protocol](protocol.md).
The full adsorption milestone is **not complete**: application packed-bed
transport fails uniform-concentration preservation. A failed coefficient is
not passed to adsorption. Stock convective-outlet controls pass separately.

## Prerequisites and exact commands

The target machine remains Apple Silicon macOS. Reuse the existing project
Python/Julia environments, pinned OpenLB 1.9r0 distribution and private DWSIM/.NET
installation. No dependency or system installation changed. Native prerequisites
are the original [synthetic bed](../README.md) N=24/32/48 fields. Existing
completed examples and their results are not overwritten.

From the repository root:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/preflight.py --output results/adsorption/column_transport/preflight.json
julia --project=julia julia/adsorption/test_column_boundary_reference.jl
julia --project=julia julia/adsorption/test_column_boundary_profile.jl
julia --project=julia julia/adsorption/test_operating_point_hydraulics.jl
julia --project=julia julia/adsorption/test_film_particle.jl
uv run python -m examples.adsorption.column_transport.hydraulics
uv run python -m examples.adsorption.column_transport.hydraulic_diagnostics
```

Use the hydraulic run directory printed by that command for the scalar study:

```sh
uv run python -m examples.adsorption.column_transport.run --hydraulic-run results/adsorption/column_transport/hydraulics/runs/20261110T154704Z-d14d13ed
uv run python -m examples.adsorption.column_transport.sensitivity
uv run python -m examples.adsorption.column_transport.process
uv run python -m examples.adsorption.column_transport.check
uv run python -m examples.adsorption.column_transport.process --check
uv run python -m examples.adsorption.column_transport.plot
uv run pytest tests/adsorption/test_column_transport.py -q
uv run python scripts/check_repository.py
```

The dated path above identifies this machine's actual hydraulic run. For a new
execution substitute its printed unique directory. Both stages create unique
run IDs and explicit `latest.json` pointers. A new run does not replace old
evidence. The scalar coordinator also accepts an already completed d800 N32
fixed-flow group in a hydraulic run still executing its independent refinements.

`hydraulic_diagnostics` normally executes two new native cases using the latest
completed hydraulic study. Its optional `--hydraulic-run PATH` selects another
run. `--reuse-native PATH` only reanalyzes a saved diagnostic run after verifying
all archived hashes and identical native binary, physical argv and input bytes;
that route executes zero new native cases. A wrapper-only portability update
was verified this way without claiming another simulator execution.

The selected final process route is `process`: actual native output is passed to
Julia's matched-Reynolds operating-point translation, and those Julia outputs
set the actual DWSIM pump inputs. `process --check` replays Julia and checks saved
native/process evidence; it does not launch DWSIM. The earlier hydraulic driver
also retains its independently checked Python arithmetic and three DWSIM cases;
the final Julia route adds two actual DWSIM executions. These are five executed
cases, not five different particle sizes.

`check` verifies saved hashes, reanalyzes native output and replays the Julia
boundary reference. It launches **zero new OpenLB or DWSIM runs**. An evidence
check passing does not convert failed physical/numerical screens into a pass.
`plot` only plots saved output and a Julia reference. No new browser/UI claim
is made. Full-project PyVista tests require host graphics access on this machine.

## Component boundaries and SI contract

`openlb/` contains the scalar application; `openlb_flow/` is a separate copy of
the established native FCC flow case, without changes to its numerical kernel.
The original vendor tree is never written. External builds stay under
`.local/adsorption/column_transport*`. All descriptors, collisions, momenta,
boundaries, geometrical indicators and data structures come from OpenLB.
Generated model-specific config copies are retained in each scalar run.

Scalar model IDs are 0: stock TRT D3Q27; 1: stock combined regularized TRT D3Q27;
2: stock MRT D3Q7; 3: stock TRT D3Q7. Regularization is the upstream combined
class, not a project collision implementation. TRT channels are unequal.
MRT uses the upstream relaxation spectrum. In the current frozen native raw
schema `tau_odd` gives its diffusive relaxation input; the algebraic `tau_even`
and `magic` fields apply to TRT models and **are unused for MRT**. Interpret
the descriptor using `model_number`, not the generic descriptor string.

Coordinates are metres, x is the positive flow direction, and lattice i/j/k
are zero-based. Node centres are `(i+.5)*dx`. Scalar ends are periodic for
constant-preservation tests; finite domains have a stock Dirichlet inlet and
either stock zero-gradient or D3Q7 Robin/Neumann exit. The finite scalar control
volume excludes reservoir planes, with length `(nx-2)*dx`. Physical pore volume
excludes solid nodes. Temporary bounce-back population storage is reported
separately and is never counted as adsorbent or physical liquid volume.

| File | Quantities |
| --- | --- |
| `request.json` | Actual argv, model, scales, native velocity and executable hashes |
| `raw.json` | dx, dt in SI, model selector, fluid volume m3, flow m3/s, timing |
| `inventory.csv` | Time s; fluid/wall numerical solute kg; cumulative inlet/outlet kg; residual kg; fluid concentration extrema kg/m3 |
| `flux.csv` | Midstep time s, signed +x inlet/exit total scalar flux kg/s |
| `final.csv` | i/j/k, material, x/y/z m, concentration kg/m3, ux/uy/uz m/s |
| hydraulic `field.csv` | Actual native physical velocity and density at fluid nodes |
| DWSIM `summary.json`/`.dwxmz` | Actual solved and save/reload checked pure-water pump |

Kinetic link flux and the same streaming control volume close the scalar mass
balance. The initial inventory is explicit, not an inlet contribution. The
stock outlet's `NEIGHBOR` field and changed boundary populations are synchronized
across periodic transverse copies using stock communicators. No clipping,
residual correction, artificial tail or geometry/inventory fitting is used.

## Guarded adsorption handoff

`bridge.adsorption_column` reuses the unchanged `FilmParticle.jl` component and
requires a traced convective-column closure whose preservation, mass, bounds,
mesh/time, outlet and longer-segment screens all passed. Its dispersion cannot
be overridden with a per-column fit. Current native evidence is rejected;
effective dispersion and new service predictions remain JSON null.

`coupling.predict` additionally freezes source equilibrium, particle-envelope
density and the batch-only fitted particle diffusivity. It enforces full dry
inventory and the native particle/voidage/velocity/diffusivity operating point.
The column film coefficient uses the declared correlation with multiplier one;
the fitted batch film coefficient is separately traceable. Hypothetical passing
packets in unit tests qualify this software contract, not an actual pore-to-column
prediction. Actual native packets remain subject to the failed-data gate.

The separate `ColumnBoundaryReference.jl` is an independent uniform 1D reference,
not another pore-scale solver. It distinguishes fixed-concentration inlet from
imposed-total-flux (Danckwerts) inlet. The Python tests also call the real frozen
Julia particle/column model in a nonadsorbing limit, checking physical void time,
mass, inventory and convergence to analytical convective-column moments.

## Hydraulic scenario and limitations

The independent hydraulic branch uses the source flow 43.8 mL/min and 0.80/0.45 mm
diameters, but synthetic FCC at nominal voidage .6. It preserves the source
density-derived 7.032455 g inventory by choosing length .14325 m at 11 mm ID and
particle-envelope density 740/(1-.427) kg/m3. This does not reinterpret printed
bulk density as particle density or claim experimental FCC geometry. DWSIM water
properties at 295.15 K set Reynolds-matched native velocities. Pump efficiency
.75 is chosen. Pressure is adjusted to impose flow, not fitted to breakthrough.

At these finite Reynolds numbers the measured quantity is operating-point
hydraulic resistance. It is not silently extrapolated as a flow-independent
Darcy permeability. This branch may verify pump duty even if solute transport
fails. A usable-capacity/pumping tradeoff still requires the missing qualified
adsorption calculation. Experimental pressure validation is unavailable.
