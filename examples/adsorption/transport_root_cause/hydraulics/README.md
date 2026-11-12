# Native fixed-flow water refinement

This component adds independently converged N64/N96 velocity fields to the
existing d800 synthetic FCC operating point. It uses the same stock OpenLB
TRT flow model and records the actual flow and voxel liquid fraction. It is
not an adsorption fit or a substitute for scalar transport convergence.

From the project root:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.hydraulics.run --grid 64
```

An optional finer run can use an explicitly named completed refinement for
its initial pressure estimate:

```sh
uv run python -m examples.adsorption.transport_root_cause.hydraulics.run --grid 96 \
  --previous-refinement results/adsorption/transport_root_cause/hydraulics/runs/RUN_ID
```

Each new resolution starts from rest. Pressure corrections on that same
resolution start from the previous converged native field, then independently
converge again. `source_point.json` is written only after all declared flow,
stationarity, density, Mach, porosity and plane mass-flux screens pass. It
contains the exact SI velocity CSV path, hashes, carrier viscosity scaling,
actual flow, interstitial speed and grid resolution for the scalar runner.

`water_population_audit.csv` observes a further stock collision/streaming step
after flow convergence. Its generic `concentration_*` columns represent native
dimensionless water density; `water_population_audit.json` records that meaning.
The audit reads populations without altering any model equations or parameters.

The [protocol](protocol.md) declares the unchanged physical conditions and
acceptance screens. Run-specific requests, inputs, binary, source snapshots,
histories and geometry records are retained under `results/`.

## Checked stock OpenMP execution

The separate `omp_check` command compiles the stock OpenLB core and this same
application with matching OpenMP flags, then fully converges the same N32 case
in serial and with four threads. On this machine the checked pair produced
identical exported velocities and used 75.37 s versus 21.95 s. This timing
comparison concerns this small case only. All native acceptance screens passed.

```sh
uv run python -m examples.adsorption.transport_root_cause.hydraulics.omp_check \
  --serial-binary results/adsorption/transport_root_cause/hydraulics/runs/20261110T205236Z-c959571b/binary/capillary3d \
  --threads 4
uv run python -m examples.adsorption.transport_root_cause.hydraulics.run_omp \
  --grid 64 \
  --crosscheck results/adsorption/transport_root_cause/hydraulics/omp_checks/20261110T205922Z-84ff1b38
uv run python -m examples.adsorption.transport_root_cause.hydraulics.run_omp \
  --grid 96 \
  --crosscheck results/adsorption/transport_root_cause/hydraulics/omp_checks/20261110T205922Z-84ff1b38 \
  --previous-refinement results/adsorption/transport_root_cause/hydraulics/omp_runs/20261110T210317Z-fdf047ba
```

These are the tested local run paths; a fresh execution prints its own new
directory. `run_omp` verifies the completed cross-check, exact binary and
hydraulic parent before using that binary. It records four threads and disables
dynamic thread counts within that process; it does not change shell or system
settings. `stock_omp.build_case_omp(source_directory, build_directory,
output_directory)` exposes the same external build approach for other diagnostic
cases, whose numerical agreement must be checked separately.

## Preserved N96 interruption and continuation

The original N96 cold solve completed and passed stationarity, mass, Mach,
density and porosity checks, but its flow was 0.536% above target. Its second
pressure trial did not finish before the Codex session update. The partial
history ends at step 38016; there is no final field or native return status for
that trial. Both original attempts and their commands remain in
`omp_runs/20261110T212155Z-83b2b0e5/`, with `interruption.json` guarding them.

The continuation reuses the first converged N96 field and exact archived
executable. The last partial velocity informs only the next pressure guess.
It supplies no accepted velocity field or stationarity observation. A new
native solve must converge and satisfy the original screens at fixed flow.
The archived velocity/density CSV is an initial estimate, not a checkpoint of
all lattice populations. The continuation does not resume the partial trial's
time index.

```sh
uv run python -m examples.adsorption.transport_root_cause.hydraulics.continue96 \
  --interrupted results/adsorption/transport_root_cause/hydraulics/omp_runs/20261110T212155Z-83b2b0e5
```

The single resumed pressure trial passed all six screens after 43,777 native
steps. The accepted reference-column flow was 7.3000006683e-7 m³/s (relative
error 9.15e-8); actual voxel porosity was 0.6002242477. The computational
carrier's pressure gradient was 2276.9603 Pa/m, 0.556% below N64. The exact
accepted field, hashes and population checks are in the
[compact N96 record](reference/grid96_reference.json). The accepted output
folder is `continuations/20261110T223616Z-3f4c960a/` under the hydraulic results.
