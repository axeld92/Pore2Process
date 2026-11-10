# Reproducing the connected hydraulic showcase

The numerical models, adapters and evidence records are included in this
repository. Full article files, installed simulators and routine generated
results are deliberately absent. The native setup below targets the verified
**macOS ARM64** configuration; other platforms have not been verified.

## Python–Julia reference first

Install uv and Julia through their official installation instructions. From
the cloned project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run python scripts/preflight.py
uv run python -m pore2process_demo.cli --output results/default
uv run pytest -q
julia --startup-file=no --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Julia must resolve to the same executable selected by bootstrap. Pass
`--julia /absolute/path/to/julia` to bootstrap if necessary. The VTK tests need
normal host graphics access. The capillary-network tab can be used immediately:

```sh
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

Open `http://127.0.0.1:5006/app`, select **Capillary network** and press
**Run simulation**. This reference does not require OpenLB or DWSIM.

## Native tool prerequisites

The tested setup reuses Homebrew GCC 15.2 and tinyxml2 11.0 at `/opt/homebrew`,
plus macOS zlib. Inspect the machine before installing or changing these
prerequisites. The current build configuration names `g++-15` and `gcc-15`
explicitly. See [the platform inventory](environment_check.md) and
[compatibility setup](../compatibility/README.md).

```sh
uv run python compatibility/bootstrap_tools.py
uv run python compatibility/bootstrap_dotnet.py
uv run python compatibility/preflight.py
```

Bootstrap downloads pinned OpenLB 1.9r0, DWSIM 10.2.15 and private .NET SDK/runtime
artifacts and verifies their hashes. Internet access and download space are
needed initially. Tool distributions and build products stay in ignored
`.local/`; original OpenLB files are audited against the archive.

## Generate the prerequisite process evidence

The newer studies validate their saved ancestors. They do not substitute the
compact documentation references for a fresh native run. From a clean checkout,
run this chain in order:

```sh
uv run python -m graded_support.run
uv run python -m graded_support.check
uv run python -m process_filtration.run
uv run python -m process_filtration.check
uv run python -m crossflow_recycle.run
uv run python -m crossflow_recycle.check
uv run python -m crossflow_benchmark.run
uv run python -m crossflow_benchmark.check
```

The graded-support stage builds the native code and obtains the DWSIM water
property anchor. Subsequent stages solve and save actual DWSIM flowsheets.
The benchmark can use tracked literature extractions without the ignored
article inbox. Add `--require-sources` only when the matching original files
are available; present article files must match their retained hashes.

These steps can take substantial time. Every attempt gets its own run folder;
`latest.json` advances only after acceptance. A failed attempt preserves its
logs and must be diagnosed before continuing. The older ancestor protocols
and numerical references remain distinct from the current TRT showcase.

## Generate the explicit-channel showcase

```sh
uv run python -m graded_support.run --protocol star_channel/support_protocol.json --output results/graded_support_trt --workers 2
uv run python -m star_channel.native_run
uv run python -m star_channel.run
uv run python -m star_channel.check
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

The explicit showcase uses freshly checked TRT support evidence, native
star-channel data and the preceding crossflow process settings. The Panel app
then offers **Filtration showcase**. Press **Run filtration process** to launch
an actual Julia/DWSIM calculation from the completed native evidence.

Do not run `star_channel.reference` merely to use the app: that maintainer
command replaces the compact tracked documentation reference. A normal new
study remains in ignored `results/` with its own snapshots and checks.

## Reproduce the closed filtration operating map

After the preceding native/channel and process chain has passed:

```sh
uv run python -m operating_map.evidence
uv run python -m operating_map.native
uv run python -m operating_map.run
uv run python -m operating_map.check
julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
uv run pytest -q
```

The native study verifies the higher-Re fully developed laminar branch; the
process study evaluates the conditional material/operation map in actual DWSIM.
This extension produces a report, figure, profiles and comparison/selection CSVs.
The Panel interface retains its previously verified low-Re controls. See
[operating_map/README.md](../operating_map/README.md) for bounds and resumption
commands, and [the final report](operating_window.md) for results.

The `operating_map.reference` command is a maintainer tool that replaces tracked
documentation references; it is not required for normal reproduction. The
filtration hydraulic case is now closed for active development and retained
for regression. The separate adsorption case is outside this publication.

## What is and is not established

The development workspace completed the full native/DWSIM sequence and browser
checks; the numerical code and final UI hashes are retained in the tracked
references. Publication verification additionally checks the Python–Julia
reference from a separate source checkout. Full native re-execution on a new
machine, other operating systems, offline browser assets and production load
are not established by publishing the files.

See [the scientific report](star_channel_showcase.md),
[acceptance history](implementation_status.md) and [browser record](browser_checks.md)
for the actual numerical ranges, conditional biological-retention scenario and
remaining research limits.
