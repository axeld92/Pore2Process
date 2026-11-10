# Python–Julia setup and run guide

This guide contains the retained setup procedure, tested runtime versions and
export contract. The maintained scientific explanation is the
[capillary-network report](../../examples/hydraulics/capillary_network/REPORT.md).
For later native/process studies, use [ordered reproduction](reproduction.md).

## Run on this Mac

From `/Users/dula/Documents/ChatGPT/Pore2Process`:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run panel serve app.py --address 127.0.0.1 --port 5006
```

Open [the local app](http://127.0.0.1:5006/app). With the verified native study
available it opens **Filtration showcase**; press **Run filtration process**.
Select **Capillary network** and press **Run simulation** for the original
Julia network/Darcy reference. Without the later study evidence, that original
tab opens by default.
Add `--show` to the serve command to request automatic browser opening.
Stop the server with Control-C. The app listens only on this computer.
The environment script adds the project-local uv installation to PATH and keeps
Python/graphics caches inside `.local/`; source it from the project root.

The selected Julia executable is:

```text
/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia
```

`/Users/dula/.juliaup/bin/julia` resolved to that same executable during setup.
Bootstrap saves its absolute path in ignored `.local/julia_runtime.json` and
sets up the single `julia/` project. Both JuliaCall override variables are set
before import, PythonCall is pinned to exactly the installed JuliaCall version,
and application launches never rewrite manifests or install Julia packages.
If the project moves, Julia changes, or the Python lock is updated, rerun
bootstrap. Use `--julia /absolute/path/to/julia` to select a different installed
Julia deliberately. The normal user Julia depot is used; no second Julia
installation or separate Julia bridge environment is created.

## Recreate setup

Install Julia using the [official installer](https://julialang.org/install/), or
reuse an existing compatible Julia (tested here: 1.11.6). Install uv through an
[official supported method](https://docs.astral.sh/uv/getting-started/installation/).
For a project-local uv installation, this machine used:

```sh
/Users/dula/.pyenv/versions/3.11.9/bin/python3 -m pip install --target .tools uv
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
```

On another machine, use an available Python 3 with pip for the first command.
uv selects Python 3.12 and creates `.venv`; it can download Python if needed.
This Mac already had Homebrew Python 3.12.9. System Python and other Python
environments were not changed. Initial package downloads require internet access.

Tested on macOS 26.6.2 / Apple Silicon with Python 3.12.9, Julia 1.11.6,
uv 0.13.7, NumPy 2.5.4, Panel 1.10.0, PyVista 0.49.2, VTK 9.7.2,
Matplotlib 3.11.2, JuliaCall/PythonCall 0.9.37, and pytest 9.2.0.
`uv.lock` and `julia/Manifest.toml` record the resolved dependencies.

## Tests and exports

```sh
source scripts/env.sh
uv run python scripts/preflight.py
uv run python scripts/preflight.py --https
julia --startup-file=no --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
uv run pytest -q
uv run python -m pore2process_demo.cli --output results/default
uv run python scripts/verify_reference.py
```

For an exact Julia selection, replace `julia` in the test command with the
absolute executable above. Preflight uses only the Python standard library,
distinguishes absent packages from failed imports, and runs bounded probes.
Its optional HTTPS metadata checks do not establish artifact availability.
It defers importing JuliaCall until bootstrap exists to avoid automatic downloads.

On this Mac, running the full UI test inside Codex's terminal sandbox aborted
when PyVista initialized macOS graphics. The full suite passed when run with
normal macOS graphics access. Numerical/CLI operations do not need an OpenGL
window. Do not mistake the sandbox graphics restriction for a solver failure.

CLI inputs use SI units, for example:

```sh
uv run python -m pore2process_demo.cli --p-in-Pa 200 --radius-m 0.00005 --output results/changed
```

Each CLI export directory and UI download contains:

| File | Contract |
| --- | --- |
| `summary.json` | Completed input snapshot, SI units, run ID/time, counts, geometry convention, flow, inferred permeability, column inputs/result, diagnostics, runtime versions, timing |
| `nodes.csv` | Zero-based node ID, x/y/z in m, pressure in Pa |
| `throats.csv` | Zero-based from/to IDs, radius and length in m, signed from-to flow in m³/s |
| `pressure_profile.png` | Numerical mean pressure and independent analytical line, labelled axes |
| `network.vtp` | Original SI line geometry with node IDs, pressure, throat flow/radius/length and orientation metadata |

Unavailable values are JSON `null`. UI ZIPs use private temporary directories
and immutable completed-run snapshots. Editing network controls does not solve
automatically; it marks old results as stale. Failed runs preserve the previous
successful result. Column edits call only the separate Julia Darcy function;
invalid column inputs disable export until corrected. Downloads retain the run
ID of the network supplying permeability and its viscosity.


## Model, architecture and recorded acceptance

See the [network report](../../examples/hydraulics/capillary_network/REPORT.md)
for assumptions, equations and results; the [architecture](../architecture/overview.md)
for common data/runtime boundaries; and the [dated acceptance record](../verification/first_simulation.md)
for the checks actually performed. The single Julia project, source snapshots
and unavailable-value conventions are retained.
