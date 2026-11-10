# Explicit hydraulic filtration in DWSIM

This milestone extends the saved fixed-opening OpenLB study into a once-through
pure-water flowsheet. Julia retains the pressure–flux law; actual DWSIM pump,
splitter, material/energy streams and pressure-loss valves solve the process.
Recovery and a separate channel pressure drop are declared scenario inputs.
No additional article or experimental measurement is required for this step.

From the project root on the tested Mac:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.filtration.once_through.run
uv run python -m examples.filtration.once_through.check
uv run pytest -q
/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
```

The runner uses the existing private DWSIM 10.2.15.0 and .NET SDK 10.0.403
installed by [the compatibility setup](../../hydraulics/tool_checks/README.md). It builds
into `.local/process_filtration/runner/`; no installed simulator files are
changed. Python and Julia use the existing `.venv`, locks and shared `julia/`
project. Julia's usual user compilation cache must be writable; a restrictive
automation sandbox needs access to that cache. No new dependencies are added.

By default the input is `results/graded_support/latest.json`. It must reference
a passed original study with its native files and actual DWSIM evidence intact.
The checker replays that upstream numerical analysis and verifies its hashes
without updating the original verification record. This step makes **18 fresh
DWSIM filtration solves** and save/reload solves, and **zero new OpenLB solves**.

To select a particular source or result:

```sh
uv run python -m examples.filtration.once_through.run --source results/graded_support/runs/20261108T222740Z-26e9f37b/summary.json
uv run python -m examples.filtration.once_through.check results/process_filtration/runs/RUN_ID/summary.json
```

Each execution creates a fresh directory and updates
`results/process_filtration/latest.json` only after success. Earlier runs are
retained. Failed runs retain status and logs. The source summary is copied
byte-for-byte; its original path, SHA-256, prior numerical checks and all
upstream provenance remain traceable. Source data are not silently regenerated.
If starting from a clean checkout without local results, first reproduce
[the upstream study](../graded_support/README.md); the compact documentation
reference is not a substitute for the native measurements.

Outputs include:

- `summary.json`, `inputs.json`, the source snapshot and numerical verification;
- `comparison.csv`, with units in field names and explicit reference volumes;
- `energy_comparison.png` and an actual DWSIM-rendered `flowsheet.png`;
- each case's actual `filtration.dwxmz`, stream summary, request, log and render;
- installed API inventory, build log, implementation/runtime hashes and timings.

The checker verifies byte integrity, Julia loads, every unit/overall balance,
liquid states, recovery, pressure levels, specific energy, saved/reloaded stream
states, original DWSIM connector IDs/ports, and cross-case scaling. It does
**not** rerun external simulators. Fresh DWSIM solves require the run command.

See [the scientific model, results and limits](../../../docs/process_filtration.md)
and [the interface decision](../../../docs/dwsim_filtration_interfaces.md).
The original Panel demonstration is independent of this command-line study.
