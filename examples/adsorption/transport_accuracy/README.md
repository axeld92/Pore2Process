# Engineering accuracy of packed-bed transport

The subsequent [bounded investigation and handoff](DIAGNOSIS.md) is complete.
It defers this high-Peclet transport branch after four additional native controls;
the previous fourteen-run assessment remains intact.

Read [REPORT.md](REPORT.md) for engineering meaning and the three prospective
protocols: [accuracy](protocol.md), [known smooth inputs](ramp_protocol.md), and
[stock regularized inlet](inlet_protocol.md). This is a separate assessment of
the existing stock model; previous strict failures remain unchanged.

## Environment and execution

Reuse the existing project Python environment, Julia project, OpenLB 1.9r0,
GCC 15.2 and saved fixed-flow hydraulic fields. No installation, dependency-lock,
vendor or earlier numerical-source change is required. The current machine is
Apple Silicon macOS 26.6.2, Python 3.12.9, Julia 1.11.6 and JuliaCall/PythonCall
0.9.37. The complete local inventory is saved under
`results/adsorption/transport_accuracy/preflight.json`.

From the repository root:

```sh
source scripts/env.sh
uv run python scripts/preflight.py --output results/adsorption/transport_accuracy/preflight.json
uv run python -m examples.adsorption.transport_accuracy.run
uv run python -m examples.adsorption.transport_accuracy.ramp
uv run python -m examples.adsorption.transport_accuracy.regularized_inlet
uv run python -m examples.adsorption.transport_accuracy.profile
```

The first command executes four native cases: two periodic uniform tests and
two sharp-inlet fronts, one per diameter. The second executes four smooth-input
cases, two input durations per diameter. The third executes two short inlet
diagnoses. The fourth executes four periodic smooth-profile controls, comparing
free flow against an exact analytical solution and porous flow against its
concentration bounds. It follows [profile_protocol.md](profile_protocol.md).
Each creates a unique directory and prints it. All use actual OpenLB.
They do not launch DWSIM, fit adsorption or modify the original distribution.
Each driver validates the saved hydraulic parent and its complete source/artifact
hashes. The first three drivers support `--hydraulic-run PATH`; the periodic
profile control uses the latest completed hydraulic parent.

For each of the first three printed native directories, execute the separate assessment:

```sh
uv run python -m examples.adsorption.transport_accuracy.assess --run PATH_TO_NATIVE_DIRECTORY
uv run python -m examples.adsorption.transport_accuracy.check
```

After all three stages are assessed:

```sh
uv run python -m examples.adsorption.transport_accuracy.check --all
uv run python -m examples.adsorption.transport_accuracy.profile --check
uv run python -m examples.adsorption.transport_accuracy.plot
uv run pytest tests/adsorption/test_transport_accuracy.py -q
uv run python scripts/check_repository.py
```

The native runs recorded on this machine are:

| Stage | Directory under `results/adsorption/transport_accuracy/` |
| --- | --- |
| Screening | `runs/20261110T180254Z-e4d0abff` |
| Known smooth inputs | `ramp/runs/20261110T180944Z-096133a1` |
| Stock inlet diagnosis | `regularized_inlet/runs/20261110T181253Z-eecee644` |
| Periodic profile controls | `profile/runs/20261110T182020Z-41c6cae7` |

All four have saved-evidence checks; the full Python suite passes 180 tests.
The [report](REPORT.md) records the small nominal-versus-actual finite-window
timing discrepancy and all failed numerical screens.

`check` verifies hashes and replays actual Julia moment analysis; it runs zero
native or DWSIM simulations. A passing saved-evidence check does not turn a failed
physical/numerical screen into a pass. Full-project graphics tests require host
graphics access on this machine. No new browser/interface acceptance is claimed.

## Files and units

Full runs stay in ignored `results/adsorption/transport_accuracy/`. Native
execution and later assessment have separate IDs. Their snapshots, commands,
build records, versions, numerical outcomes and source hashes are retained.
Curated [reference records](reference/README.md) link to exact summaries and hashes.

The baseline reuses the unchanged `column_transport/openlb/` application, built
separately. `openlb_ramp/` changes only the prescribed input schedule, using the
stock `PolynomialStartScale`. `openlb_regularized_inlet/` changes only that
inlet's stock boundary selection and adds extrema-location diagnostics.
`openlb_profile/` changes only the analytical initial profile and requires
periodic ends; it implements no transport algorithm. Stock
TRT, regularization, descriptor, geometry, walls, data structures and transport
operators remain upstream OpenLB facilities. External builds are under
`.local/adsorption/transport_accuracy*`; none is placed in the vendor tree.

Native coordinates and zero-based indices retain the earlier contract. X is
flow direction; positions are metres, velocity m/s, concentrations kg/m3, water
flow m3/s, solute flux kg/s, inventories kg, diffusivity m2/s and time seconds.
The unit-cell water Q is its own cross-section times superficial velocity; it is
scaled from the same column superficial duty, not the whole-column 43.8 mL/min
through a microscopic cell. Finite physical liquid volume excludes solid nodes
and reservoir planes. Numerical wall-population storage is reported separately.

`feed.csv` records the actual prescribed ramp at each time-interval midpoint;
`flux.csv` gives unnormalized inlet/outlet solute fluxes. `inventory.csv` contains
physical-liquid inventory, numerical wall storage, cumulative actual fluxes and
liquid extrema. The regularized-inlet variant adds the i/j/k locations of those
extrema. `final.csv` retains full fields. All JSON unavailable metrics are null.

For a known ramp input the response moments include its delay/spreading. Julia
quadrature evaluates the recorded input and output separately; input mean and
variance are subtracted using the convolution identity, without clipping,
deconvolution, fitted delay or invented tail. The stock quintic ramp's independent
input moments are T/2 and T²/28. Correcting the input does not qualify an unphysical
response. No failed coefficient reaches the historical guarded particle bridge.

## Later stages

`run --stage refine --parent PATH` supports N32 half-dt and N48 at that same
physical dt. `run --stage holdout --parent PATH` supports a doubled segment.
The ramp drivers support `--n 48 --half-dt`, `--cells 2`, and longer windows.
These are documented diagnostic capabilities, not claims that those stages
passed. The protocols require bounded transport before accepting a closure and
the engineering service-error budget. Earlier batch and pumping results retain
their separate validity and scope.

## Bounded diagnosis, now completed

This follow-up reuses the exact preserved profile executable and changes only
declared diagnostic inputs. It compares stationary water, half/tenth dt and a
matched-dt N48 field. It does not rebuild OpenLB, fit an adsorption curve or
launch DWSIM. See [diagnosis_protocol.md](diagnosis_protocol.md) for its stopping
rule and [DIAGNOSIS.md](DIAGNOSIS.md) for the conclusion.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_accuracy.diagnosis
uv run python -m examples.adsorption.transport_accuracy.diagnosis --check
uv run python -m examples.adsorption.transport_accuracy.plot_diagnosis
uv run pytest tests/adsorption/test_transport_diagnosis.py tests/adsorption/test_transport_accuracy.py -q
```

The recorded native run is
`results/adsorption/transport_accuracy/diagnosis/runs/20261110T183934Z-2a431270/`.
The checker validates 27 source files and 573 native artifacts and replays the
field audits, with zero new simulator executions. Ten targeted tests pass.
