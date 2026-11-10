# Synthetic fixed-flux filtration study

[Scientific report: purpose, theory, results and conclusions](REPORT.md).

Read [the research scope and contracts](../../../docs/verification/filtration/layered_support.md).
This workflow adds a fixed selective layer to fresh OpenLB support simulations,
then calculates an actual DWSIM pumping duty. It is not experimentally calibrated.

From the project root on the configured Mac:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.filtration.layered_support.run
uv run python -m examples.filtration.layered_support.check
uv run python -m examples.filtration.layered_support.validate_inputs
```

Setup reuses the existing Julia 1.11.6/PythonCall project, Python lock, OpenLB
1.9r0 distribution, GCC 15.2 and local DWSIM 10.2.15/.NET installation. No new
dependency is needed. `filtration.build` uses the shared external build helper;
the native case includes `examples/hydraulics/capillary_benchmark/openlb/case.h` and a small application
driver that changes dimensions/input bounds. Native build products go to
`.local/filtration/build/`. All original OpenLB sources are checked against the
pinned archive before and after building. The original capillary case retains
its defaults and commands.

The workflow creates a fresh `results/filtration/runs/<run-id>/` and updates
`results/filtration/latest.json` only after acceptance. Failed partial runs are
retained with their logs; they are not marked passed. Outputs include:

- Per-support/per-grid actual OpenLB `raw.json`, `planes.csv`, `midplane.csv`,
  `history.csv`, analysis and logs; low/doubled-pressure and doubled-domain cases.
- Every DWSIM input/result and saved `water_pump.dwxmz` flowsheet.
- `summary.json`, `comparison.csv`, `filtration_comparison.png` and provenance
  hashes, numerical settings, units, versions, cold Julia initialization and
  actual solver/adapter timings.
- A separate `verification.json` after the check command recomputes the Julia
  analysis and validates stored exports. This check reads DWSIM evidence; it
  does not rerun DWSIM.

To reuse verified native measurements but perform fresh process solves:

```sh
uv run python -m examples.filtration.layered_support.run --reuse results/filtration/runs/<run-id>/summary.json
```

Reuse requires verified file hashes, the same solver binary, study settings
and fluid properties. Normal execution always performs fresh native solves.
Keep the full run directory when moving a snapshot, rather than only its summary.
The supported operating settings are currently explicit constants in
`examples/filtration/layered_support/run.py`; there is no new UI or arbitrary-geometry input contract.

Julia's independent numerical tests and Python integration tests:

```sh
julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
uv run pytest -q
```

The first command must resolve to the same Julia as the bridge. On this Mac,
use `/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia`
explicitly if needed. The Codex sandbox needs normal user-cache access to
precompile Julia and normal macOS graphics access for the original PyVista UI
tests; this is separate from the numerical model. No new package installation,
browser check or connected UI is implied by this CLI milestone.

Article Markdown/PDFs and supplements can be supplied locally in
`literature/inbox/`. See [the exact article request](../../../literature/README.md).
