# Controlled graded-support study

[Scientific report: purpose, theory, results and conclusions](REPORT.md).

[Verified results and figures](../../../docs/verification/filtration/graded_support.md): the widest
interior reduces native resistance by 35.88% and conditional transmembrane
shaft energy by 11.30%. The final saved reference is
`results/graded_support/runs/20261108T222740Z-26e9f37b/summary.json`.

Three synthetic connected pores use the same 4-µm opening radius, 40-µm body
length, 20-µm representative cell width and fixed terminal necks. Interior
radii are 4, 5 and 6 µm. The [protocol](protocol.json) defines the geometry,
checks and literature scaling before the simulations. This extends the
completed examples in a separate folder and leaves the original OpenLB tree
untouched.

## Run on the configured Mac

From `/Users/dula/Documents/ChatGPT/Pore2Process`:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.filtration.graded_support.run
uv run python -m examples.filtration.graded_support.check
uv run python -m examples.filtration.graded_support.validate_inputs
uv run pytest -q
/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
```

Setup reuses CPython 3.12.9, Julia 1.11.6, matching JuliaCall/PythonCall 0.9.37,
OpenLB 1.9r0, GCC 15.2 and the installed DWSIM/.NET adapter. No new dependencies
are required. Native cases compile into `.local/graded_support/build/` using
the existing external build helper and stock core. Its vendor audit checks
2,070 original files against the pinned distribution archive.

Two independent native subprocesses run at a time. `--workers 1` runs them
serially. Worker threads only launch subprocesses in private directories;
Julia calls and analysis stay on the calling main thread. This is a CLI study,
not a background Julia worker or a new interface.
Each native process has a 30-minute execution allowance; the shared launcher's
default for earlier examples remains ten minutes. Numerical iteration limits
and all frozen acceptance settings are unchanged.

Each invocation creates an immutable `results/graded_support/runs/<run-id>/`.
`latest.json` is written only after acceptance and replay. To reuse intact
native measurements from a previous or partially completed run while
repeating analysis and actual DWSIM solves:

```sh
uv run python -m examples.filtration.graded_support.run --reuse-native results/graded_support/runs/RUN_ID
uv run python -m examples.filtration.graded_support.check results/graded_support/runs/RUN_ID/summary.json
```

The cache requires identical native binary, inputs and all measurement hashes.
Missing native cases are solved afresh. Completed outputs are never overwritten.
Failed/partial runs keep their logs and are not published as the latest pass.
Normal execution performs fresh OpenLB and DWSIM solves. The separate checker
recomputes stored native analysis and Julia loads and checks DWSIM evidence;
it does not rerun either external solver.

## Scientific and software boundaries

1. **Geometry:** the body has axial knots 0, 8, 16, 24, 32 and 40 µm and radii
   4, 4, bulk, bulk, 4, 4 µm. OpenLB's existing cylinders and cones are combined
   by its stock union operator. Numerical development sections are 16 µm at
   each end and are excluded from the body pressure-tap span. The domain
   check extends them to 24 µm without changing the support body.
2. **OpenLB:** D3Q19, stock TRT (magic parameter 3/16), Bouzidi walls, local pressure boundaries and
   ValueTracer convergence are stock facilities. The case follows mesh,
   geometry, lattice, initialization, simulation and result stages. Three
   grids resolve the fixed neck diameter with 24, 32 and 40 nodes. Every design
   also has 5/20/−10-Pa, extended-lead and half-time-step checks at N24.
3. **Julia:** measured body pressure drop and signed flow determine body
   resistance and axial permeability. A separate lubrication integral is an
   approximate reference for tapered sections and an exact straight-pore
   oracle. It never supplies the OpenLB solution.
4. **Scale translation:** one common calibration anchors the measured straight
   baseline to Zhu's reported 3.12 m/h support flux at 1 bar. Each measured
   resistance ratio changes that support impedance. Zhu PDF-v2 effective
   non-support impedance stays fixed across geometries. This assumes that the
   resolved shape effect transfers through the same unresolved-structure
   multiplier; it is not an absolute prediction for a reconstructed ceramic.
5. **DWSIM:** its liquid-water properties configure native simulations. Julia's
   load then sets a real pump discharge pressure at fixed 100 LMH, 1 m² and
   75% efficiency. Nine actual flowsheets are solved, saved, reloaded and
   solved again, covering three geometries and three fixed-impedance readings.

Unreported experimental temperature/viscosity are not invented. The DWSIM
temperature is a declared 298.15-K computational scenario. Keeping the fitted
experimental impedance in that scenario is an explicit approximation; this
does not perform a temperature correction. Reading sensitivity is not
experimental uncertainty. The [project verification scope](../../../docs/architecture/verification_scope.md)
requires no new experiments or inaccessible raw datasets.

The native geometry's nominal-flux Reynolds number is reported separately.
It is not a pore Reynolds number for the unknown experimental structure behind
the baseline calibration. The native pressure-range checks establish creeping
flow for the resolved geometry; linear response of the literature-anchored
scenario remains a model assumption.

## Exports and contracts

- `inputs.json` preserves the full protocol, literature reference, fluid
  anchor and code/lock/binary hashes. `summary.json` contains all native and
  process results, timing and acceptance checks.
- Per native case: `raw.json`, `planes.csv`, `midplane.csv`, `history.csv`,
  `stdout.log`, `measurement.json` and `analysis.json`. SI field suffixes are
  explicit. Positive flow is +x; OpenLB pressure is gauge Pa.
- Coordinates use zero-based lattice indices: x=i·dx, y=(j+½)·dx, z=(k+½)·dx.
  Each native mesh spans 12 µm transversely and is centred at (6,6) µm; the
  declared representative area is the separate 20 × 20 µm cell convention.
  The finest-grid midplane has a design-dependent number of rows, one per
  fluid node, with int64 indices and float64 SI fields.
- Fine-grid `geometry/voxels.csv` exports actual OpenLB material labels for
  six-neighbour connectivity, full radius-profile and inlet/outlet audits. No Python flow solver
  is used.
- `comparison.csv` contains nine labelled process rows. DWSIM inputs use
  absolute discharge pressure = inlet absolute pressure + Julia's load and
  mass flow = shared inlet density × filtrate volume flow. DWSIM's internal
  kW/kJ outputs are converted to W/J by the existing adapter.
- `geometries.png`, `comparison.png` and three SI `*_wall.vtp` files provide
  inspection views. VTP x omits numerical leads; y/z are shifted by +4 µm to
  the centre of the 20-µm representative cell. Colour is interpolated section
  mean pressure, explicitly not a wall stress or wall-pressure solution.
- Process subdirectories retain input/result JSON, logs and actual saved
  `water_pump.dwxmz` flowsheets. `verification.json` distinguishes replay from
  new external solves. Strict JSON contains no NaN/Infinity.

The energy metric is transmembrane pressure-load shaft energy per filtrate
volume. Crossflow circulation, fouling, backwash, osmotic pressure, motor
losses, selectivity and mechanical certification are excluded. An explicit
filtration equipment/separation/recycle flowsheet remains a later process
extension; the current DWSIM cases contain actual pumps and streams.

## Primary implementation references

- Installed `examples/laminar/poiseuille3d/case.h` for the OpenLB lifecycle,
  stock dynamics, pressure and Bouzidi boundaries.
- Installed `examples/turbulence/venturi3d/` and
  [OpenLB geometry documentation](https://openlb.gitlab.io/release/content/03_geometry.html)
  for combining stock cylinder/cone indicators.
- Existing [DWSIM adapter](../../../src/pore2process_demo/runtime/dwsim_runner/Program.cs) and
  [Julia bridge](../../../src/pore2process_demo/bridge.py), reused without modification.
