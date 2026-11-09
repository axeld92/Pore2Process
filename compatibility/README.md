# OpenLB / DWSIM compatibility check

Verified on 8 November 2026 on this Apple Silicon Mac. OpenLB 1.9r0 builds
and runs its analytical channel benchmark. Python invokes a private .NET
adapter that creates, solves, saves, reloads, and solves a real DWSIM 10.2.15
water-pump flowsheet. The original Julia/Panel example is unchanged.

These scripts remain independent compatibility probes. The separate
[connected workflow](../coupling/README.md) passes OpenLB-derived properties
through Julia to DWSIM. See [the probe results and limitations](../docs/compatibility_check.md)
and [the integration milestones](../docs/integration_roadmap.md).

From `/Users/dula/Documents/ChatGPT/Pore2Process`:

```sh
source scripts/env.sh
uv run python compatibility/bootstrap_tools.py
uv run python compatibility/bootstrap_dotnet.py
uv run python compatibility/build_openlb.py
uv run python compatibility/run_openlb.py
uv run python compatibility/run_dwsim.py
uv run python compatibility/preflight.py
```

The Python scripts use only the standard library; `.venv/bin/python` can replace
`uv run python`. Bootstrap verifies pinned archive hashes and reuses cached
distributions. Initial downloads need internet access. Archives, extracted tools,
SDK, .NET package/build caches, and compiler outputs live in ignored `.local/`.
Logs, input snapshots, JSON results, and saved flowsheets live in ignored
`results/compatibility/`. The small measured reference is committed in
`docs/compatibility_reference.json`.

This configuration targets macOS ARM64. It reuses existing Homebrew GCC 15.2.0
and tinyxml2 11.0.0, plus the macOS zlib. If these are absent, inspect the machine
and install supported compiler/dependencies before building. On this Mac their
paths are `/opt/homebrew/bin/g++-15`, `/opt/homebrew/bin/gcc-15`, and
`/opt/homebrew/include/tinyxml2.h`. The build is serial `CPU_SISD`, double
precision, C++20. It does not use the discovered MPI installation. The current builder compiles byte-identical copies of the stock example in
`.local/compatibility/openlb_example_build/`, using the shared external core
archive from `.local/hydraulic/build/core/`. Configuration and build outputs
stay outside the original OpenLB distribution. The vendor `config.mk` configured
by the earlier probe is ignored and never rewritten. Original files are checked
against the pinned archive. See [the code policy](../docs/openlb_code_policy.md).

The .NET SDK is pinned to 10.0.403 / runtime 10.0.14. The adapter has no NuGet
package dependencies: it references DWSIM assemblies in the verified portable
distribution. Python launches it as a separate process with a 180-second timeout.
The adapter accepts absolute distribution, input JSON, and output directory
arguments. JSON field suffixes declare SI units; efficiency is a fraction at the
boundary and is converted to percent for DWSIM. DWSIM kW and kJ/kg outputs are
converted to W and J/kg. Assembly resolution and the working directory point at
the portable distribution. The DWSIM desktop GUI is not launched.

`run_openlb.py` requires both inlet/outlet material counts to equal N−1,
finite analytical errors, reported convergence, velocity relative L2 error <1%,
pressure relative L2 error <5%, and lower errors after N=32→64 refinement.
`run_dwsim.py` checks liquid state, outlet pressure, mass/energy balances, pump
power against ṁΔp/(ρη), energy-stream agreement, input snapshots, saved/reloaded
flowsheets, and a factor-two power change when pressure rise doubles.

To use the adapter independently after the build:

```sh
export DOTNET_ROOT="$PWD/.local/compatibility/dotnet"
export DOTNET_CLI_HOME="$PWD/.local/compatibility/dotnet-home"
export DOTNET_CLI_TELEMETRY_OPTOUT=1
export DOTNET_GENERATE_ASPNET_CERTIFICATE=false
export DYLD_LIBRARY_PATH="$PWD/.local/compatibility/dwsim"
"$DOTNET_ROOT/dotnet" \
  .local/compatibility/runner/bin/DwsimProbe/debug/DwsimProbe.dll \
  "$PWD/.local/compatibility/dwsim" \
  "$PWD/results/compatibility/dwsim/baseline/inputs.json" \
  "$PWD/results/compatibility/dwsim/manual"
```

Prefer the Python runner, which also sets DWSIM's working directory, builds the
adapter, captures logs, and validates results. No global .NET or Mono installation,
COM registration, privileged system installation, or remote service is needed.
