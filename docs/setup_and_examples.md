# Setup and examples

A local vertical slice: Python configures a synthetic capillary lattice, Julia
solves the sparse pressure equations, and Panel displays the network and a Darcy
column demonstration. The numerical package also works directly from Julia.
No external dataset or Python production solver is used.

The new [explicit-channel and biological-feed showcase](star_channel_showcase.md)
adds a **Filtration showcase** tab to the same Panel application. Its
[ordered launch commands](../star_channel/README.md) build separate native TRT
support/channel evidence before running actual DWSIM recycle cases. The
biological case reports clean initial hydraulics and an illustrative published
COD-retention sensitivity; sustainable flux at the new geometry/speed is
unavailable. Older numerical results below remain preserved as historical studies.

The final workflow requires OpenLB, Julia, Python, and DWSIM. The first example
below establishes the Python–Julia reference. Native OpenLB/DWSIM compatibility
checks have since passed on this Mac; see
[compatibility commands](../compatibility/README.md),
[verified results and limitations](compatibility_check.md), and
[integration milestones](integration_roadmap.md).

The next hydraulic milestone is also implemented: a real 3D OpenLB capillary
solve with Python exports and Julia permeability extraction. See
[hydraulic launch commands](../hydraulic/README.md) and
[the verified benchmark](hydraulic_benchmark.md). This remains a separate
CLI workflow from the first Panel interface.

The connected hydraulic comparison runs **OpenLB → Python/Julia → DWSIM**:
two pore radii, measured permeability, a Julia Darcy pressure requirement, and
actual saved DWSIM pump flowsheets. See [launch and verification commands](../coupling/README.md)
and [the OpenLB implementation policy](openlb_code_policy.md).

Candidate research directions, including European industrial showcases and
reacting/non-isothermal bioprocess options relevant to DTU PROSYS, are collected
in [research showcase options](research_showcase_options.md).

An experimental support-permeability calibration and held-out check is now
implemented in Julia/Python. It improves prediction error but **fails the
declared 25% screen**; it does not validate an OpenLB material geometry or the
complete membrane. See [results and limitations](experimental_calibration.md)
and [run/replay commands](../calibration/README.md).

The supplied Zhu (2021) matched support/membrane series has since been fitted
and tested through Julia. The added original PDF improves the Figure 6 reading:
nominal held-out maximum error is 4.56%, and the conservative reading-corner
maximum is 7.88%, both within the unchanged 10% screen. See [the PDF refinement](zhu_pdf_refinement.md);
the earlier JPEG reading and its failed sensitivity screen remain preserved.
This is a bounded effective-hydraulic comparison using published evidence.

This is a computational-tool project. Synthetic or literature-informed
geometries, numerical verification and published-data comparisons are the
scope. New experiments, unavailable raw data and reconstruction of a particular
specimen's 3D pore structure are not prerequisites. See [verification scope](verification_scope.md).

The selected research question is now **support geometry versus pumping energy
at fixed filtration flux and unchanged selective-layer resistance**. The first
synthetic application study adds a Julia layer-resistance model and real DWSIM
pump cases to fresh OpenLB support solves. See [scope and assumptions](filtration_research_case.md),
[launch commands](../filtration/README.md), and [the literature evidence](../literature/README.md).

The controlled **fixed-opening interior-widening study** is now complete:
three geometries, 24 verified native OpenLB cases, Julia scale translation and
nine actual saved DWSIM pump cases. Widening the interior radius from 4 to
6 µm reduces native support resistance by **35.88%** and conditional
transmembrane shaft energy by **11.30%**. See [results, geometry views and limits](graded_support_study.md)
and [fresh-run/replay commands](../graded_support/README.md). The literature
baseline calibration is common to all designs; this remains a synthetic
geometry study.

The **explicit hydraulic filtration flowsheet** is also complete: stock DWSIM
pump, splitter and pressure-loss valves, with feed/permeate/retentate streams.
All 18 geometry/recovery/channel-loss cases pass balances and saved/reloaded
solves. At 50% recovery, an assumed 20 kPa upstream loss reduces the widest
support's shaft-energy saving to **3.81%**. See [results and assumptions](process_filtration.md),
[interface investigation](dwsim_filtration_interfaces.md), and
[run/check commands](../process_filtration/README.md). This is a once-through
pure-water hydraulic surrogate.

The **crossflow recycle and circulation-energy extension** is now verified:
25 actual DWSIM recycle cases, separate feed/circulation pumps, temperature
control, purge, balances and saved/reloaded flowsheets. At 50% recovery and
module flow five times permeate flow, circulation is **83.09%** of shaft duty
and the widest support saves **1.91%** of total shaft energy. Higher circulation
dilutes the material benefit further under the declared channel-resistance law.
See [results and scope](crossflow_recycle.md) and
[reproduction commands](../crossflow_recycle/README.md). Cooling heat and shaft
work remain separate.

The **literature-constrained crossflow benchmark** now adds that performance
role: real Julia flux constraints and 32 actual DWSIM cases, with explicit
source-flow conventions and passive-solute balances. At fixed 100 LMH, the
table policy selects 0.9 m/s among tested points; a conservative source-reading
policy selects 1.4 m/s. The widest support saves **7.11–10.99%** at these
conditional selected points. See [results and evidence limits](crossflow_benchmark.md)
and [run/check commands](../crossflow_benchmark/README.md). No absolute module
calibration or physical fouling-time prediction is inferred.

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

## Model and independent analytical reference

Nodes are junctions without storage or additional pressure loss. Each edge is
a cylindrical capillary of length `a` and radius `r`, with Newtonian viscosity
`mu`. Flow is steady, incompressible, fully developed and laminar, without
gravity, entrance losses, or periodic boundaries. Transverse boundaries are
closed because no edges cross them.

For logical zero-based `(i,j,k)`, node ID is `i + nx*(j + ny*k)` (x fastest).
Coordinates are `(i*a, (j+0.5)*a, (k+0.5)*a)`. In node order, existing positive
x, y and z edges are emitted once, in that order. Positive flow goes from the
first endpoint to the second. Julia uses one-based endpoint IDs internally;
the bridge subtracts one exactly once and copies arrays into ordinary NumPy
`float64`/`int64` arrays. Returned snapshots do not share mutable simulation state.

The solver assembles the sparse conductance Laplacian from
`g = pi*r^4/(8*mu*a)` and solves `K_II*p_I = -K_IB*p_B` with inlet/outlet
Dirichlet pressures. Signed throat flow is `q = g*(p_from-p_to)`.
Nodal net outflow adds `q` at the first endpoint and subtracts it at the second.
`Q_in` sums inlet nodal net outflows; `Q_out` negates their outlet sum, so
transverse boundary edges cancel without double counting. `nx=2` has no
interior solve. Zero pressure drop is assigned the exact constant state.

Independently of the sparse assembly: a path along x contains `nx-1`
capillaries in series, so its hydraulic resistance is `(nx-1)/g`.
There are `ny*nz` identical paths in parallel. The linear pressure field has
equal transverse pressures, satisfying every transverse edge and interior
balance. Consequently:

```text
L = (nx - 1)*a
A = ny*nz*a^2
p(x) = p_in - (p_in - p_out)*x/L
Q_exact = ny*nz*pi*r^4*(p_in - p_out)/(8*mu*L)
k_exact = mu*L*Q_exact/(A*(p_in-p_out)) = pi*r^4/(8*a^2)
```

Each transverse node represents a cell of area `a²`. The y/z domains are
`[0,ny*a]` and `[0,nz*a]`, leaving a half-cell beyond outer node centres.
This convention also supports a single-capillary chain.
The analytical expressions are used only for tests and comparison plots;
the production answer comes from Julia's sparse equations.

Inferred permeability is `mu*L*Q_in/(A*delta_p)`. It remains positive on
pressure reversal and is unavailable at zero pressure drop. The independent
Julia column function uses `Q_column = k*A_column*delta_p_column/(mu*L_column)`.
The column assumes the same repeating medium at a larger scale; this is an
upscaling demonstration, without empirical validation for a real material.

Default inputs are `(6,4,4)`, `a=1e-3 m`, `r=50e-6 m`, `mu=1e-3 Pa s`,
inlet/outlet pressures 100/0 Pa. Expected flow is 47.1238898038469 µL/min,
permeability 2.4543692606170265e-12 m², and column flow 0.736310778185108 mL/min
for `L_column=.02 m`, `A_column=1e-4 m²`, `delta_p_column=1000 Pa`.
See [verified numerical observations](verified_reference.json).

## Scope and implementation notes

UI dimensions are bounded to `2 <= nx <= 20`, `1 <= ny,nz <= 8`. Inputs must
be finite, dimensions integers, spacing/radius/viscosity positive, and diameter
smaller than spacing. The 3D view uses physical tube radii (1×); coordinates
are converted to mm only for display. Export geometry stays in metres.
The pressure colorbar is labelled in Pa.

One Julia runtime is shared by the local server process; controls and result
snapshots belong to each browser session. Julia runs synchronously on the main
server thread. Panel/Bokeh event callbacks use async wrappers to retain event-loop
execution; there is no `await` during a Julia solve, no background Julia worker,
and no cancellation. Bokeh 3.10 moves synchronous callbacks onto workers, so
plain synchronous callbacks are unsuitable for this bridge. A bridge guard
rejects worker-thread calls before entering Julia.

Initialization and first compilation take seconds; warmed bounded solves were
milliseconds or less in the recorded run. Plots/serialization/export add their
own time. This is not a speedup or performance claim. Timing fields distinguish
runtime initialization, bridge calls including first compilation/conversion,
core function time, and reduced linear solve time.

If 3D creation raises a Python exception, the app keeps the numerical result
and pressure plot and reports the graphics failure. Interactive 3D must still
be verified on each target platform. This prototype excludes images, extracted
networks, heterogeneous radii, transients, optimization, authentication, cloud
hosting, and production job orchestration.

See [environment findings](environment_check.md) and
[acceptance status](implementation_status.md) for the exact checks and
remaining limitations. The supplied handoff is preserved unchanged.

## Primary documentation

- [JuliaCall existing environment requirements](https://juliapy.github.io/PythonCall.jl/stable/juliacall/#Using-existing-environments)
- [Julia project environments](https://pkgdocs.julialang.org/v1/environments/)
- [Julia sparse arrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/)
- [uv projects and lockfiles](https://docs.astral.sh/uv/guides/projects/)
- [Panel VTK pane](https://panel.holoviz.org/reference/panes/VTK.html) (the linked page returned 404 during this work; installed Panel's VTK docstrings/source and published reference examples were also inspected)
- [Panel FileDownload](https://panel.holoviz.org/en/docs/latest/reference/widgets/FileDownload.html)
- [PyVista](https://docs.pyvista.org/)
