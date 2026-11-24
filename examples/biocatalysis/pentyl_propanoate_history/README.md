# Local catalyst exposure history in the packed bed

[Engineering report](REPORT.md) · [Frozen protocol](protocol.json)

The stationary catalyst retains its own cumulative exposure at each axial
location. The liquid progress is advected and dispersed; catalyst damage is
not. The original batch rate, damage and recovery parameters, 6 g inventory,
bed dimensions and film assumptions are unchanged. This is a new transient
extension of the earlier steady finite-volume calculation, not a fitted
reactor kinetic law or a new OpenLB simulation.

## Run and verify

From the repository root with the existing [isolated environments](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_history/test_model.jl
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_history/test_history.py
```

The prediction campaign contains 14 prespecified model/scenario combinations,
40 reactor conditions each. A one-process run produces one partition:

```sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.run \
  --output results/biocatalysis/pentyl_propanoate_history/new_run/part_0
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.aggregate \
  --partitions results/biocatalysis/pentyl_propanoate_history/new_run/part_0 \
  --output results/biocatalysis/pentyl_propanoate_history/new_run/combined
```

For the verified parallel run, six separate Python processes each initialize
their own Julia runtime. Tasks are partitioned by deterministic case index;
each process owns its results directory. Do not share one JuliaCall runtime
between background threads. For example, in the project root on macOS:

```sh
for part in 0 1 2 3 4 5; do
  uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.run \
    --output "results/biocatalysis/pentyl_propanoate_history/new_parallel/part_${part}" \
    --part "$part" --parts 6 &
done
wait
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.aggregate \
  --partitions results/biocatalysis/pentyl_propanoate_history/new_parallel/part_0 \
               results/biocatalysis/pentyl_propanoate_history/new_parallel/part_1 \
               results/biocatalysis/pentyl_propanoate_history/new_parallel/part_2 \
               results/biocatalysis/pentyl_propanoate_history/new_parallel/part_3 \
               results/biocatalysis/pentyl_propanoate_history/new_parallel/part_4 \
               results/biocatalysis/pentyl_propanoate_history/new_parallel/part_5 \
  --output results/biocatalysis/pentyl_propanoate_history/new_parallel/combined
```

Use exactly one of these workflows, with fresh output paths; the runners do
not overwrite previous experiment records. The aggregate command verifies
all partition and profile hashes and requires all 560 qualified histories.
Only after the joint prediction ledgers have been saved does it join the
40 experimental outlet conversions. It does not infer a catalyst age or
choose a kinetic law from that comparison.

Cold-repeat the common refinement decisions and rebuild compact evidence
(here for the sequential path; substitute `new_parallel` for the parallel path):

```sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.replay \
  --combined results/biocatalysis/pentyl_propanoate_history/new_run/combined \
  --output results/biocatalysis/pentyl_propanoate_history/new_run/replay.json
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.retain \
  --combined results/biocatalysis/pentyl_propanoate_history/new_run/combined \
  --replay results/biocatalysis/pentyl_propanoate_history/new_run/replay.json \
  --output results/biocatalysis/pentyl_propanoate_history/new_run/compact
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.plot \
  --reference results/biocatalysis/pentyl_propanoate_history/new_run/compact \
  --output results/biocatalysis/pentyl_propanoate_history/new_run/figures
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_history.check
```

## Data contract

`SpatialCatalystHistory.simulate` in [model.jl](model.jl) takes the frozen
`CompositionKinetics.Parameters`, temperature in °C, nominal feed ratio,
dimensionless apparent equilibrium constant and sample times in min.
Keyword inputs are catalyst mass kg, interbead liquid L, volumetric flow
L/min, acid and alcohol concentrations mol/L, film capacity mol/(L min) per
unit progress difference, `eta_internal` and `peclet` dimensionless, number
of cell centers and a maximum time step in min. It returns Julia arrays and
independent species-balance diagnostics; it is usable without Python.

Each compressed `profiles/<case>.npz` in the untracked partition results
contains:

- `times_min` and `cell_center_fraction`: axes in min and z/L, starting at
  inlet on the left. Cell index starts at zero in exported NumPy data.
- `bulk_progress`, `surface_progress`, `cumulative_hazard`, `activity`,
  `local_rate_mol_L_min`, `local_hazard_min_inv`: arrays of shape
  `(axial_cells, sample_times)`. Progress and activity are fractions; H is
  integrated dimensionless damage at a fixed cell; rate and hazard units
  are in their names. These describe axial averages, not bead-scale maps.
- `outlet_conversion`: dimensionless fraction at each sample time. Four
  other vectors give `generated_mol`, `effluent_product_mol`,
  `backdiffused_product_mol`, `liquid_product_mol`. All are ester equivalents
  for the single balanced reaction. Initial product inventory can be nonzero
  only in the explicitly hypothetical equilibrated-prefill case; it is not
  counted as newly generated ester.

The [reference summary](reference/summary.json) retains complete-condition
scores and numerical limits; `history.csv` preserves all sampled histories,
`bed_comparison.csv` preserves all held-out predictions and actual measured
conversion, and `numerics.csv` records the accepted mesh and time-step
refinement per case. `selected_axial_profiles.csv` uses conservative
32-bin averages of selected saved full fields. `reaction_zone.csv` gives
the reaction-weighted axial location only when its instantaneous production
is non-negligible. Reading a curve at an unknown experimental age is not
the same as measuring a transient profile. The time-mean predictions are
an explicitly prespecified summary scenario, not a fitted measurement time.

Tested with Python 3.12.9 and Julia 1.11.6 on Apple Silicon macOS. The six
prediction groups completed successfully; the longest took about 392 s
after runtime loading. The independent cold replay repeated all refinement
decisions and all saved numerical arrays exactly on 36 predeclared cases
(every model/scenario in both feeds at 80 °C and 3 mL/min, plus eight
representative moderate/slow cases). It took about 115 s after setup.
All 104 direct Julia assertions, five new Python tests and seven preceding
composition/equimolar integration tests passed. All previous source records
remain unchanged. No new external package or experimental data was required.
The code runs in project space and doesn't modify OpenLB or a saved DWSIM case.
