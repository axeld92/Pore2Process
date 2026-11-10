# Pore2Process: first Python–Julia simulation

Prepared: 7 November 2026. Status: implementation plan and verified numerical reference; the hybrid application has not been built.

## Start the new local Codex chat with this

Place this file in the intended local project directory, open that directory in Codex, and use:

> Read `PORE2PROCESS_FIRST_SIMULATION_CODEX_HANDOFF.md` in full and execute its plan. First inspect the local workspace, applicable AGENTS.md instructions, operating system, runtimes, and installed dependencies. The environment observations in the handoff came from a separate hosted workspace; verify this machine independently. Then set up isolated project environments, implement and test the Julia solver, connect it to Python, and build the small local Panel interface. Continue through the acceptance checks and document the exact launch commands and any genuine blockers. Use available terminal, documentation, testing, and browser tools where they help. Make routine implementation choices yourself within the scope below, preserve unrelated files, and report only work actually verified.

## 1. Purpose and decisions already made

Pore2Process will explore the connection between pore-scale behavior and process-scale performance. The user is a chemical engineer, computational scientist, and PhD student who likes Julia. We discussed Julia-only, Python-only, and hybrid approaches and chose to investigate **Python for the interface and integrations, with Julia for custom numerical models**.

This first example should prove the complete local path:

1. Configure a synthetic pore network in Python.
2. Solve its steady pressure-driven flow in Julia.
3. Return pressures, throat flows, and effective permeability to Python.
4. Inspect the results in a local browser interface.
5. Use that permeability in a simple Darcy column calculation and export reproducible results.

Use **Python + Panel + JuliaCall + Julia**, with **PyVista/VTK** for the network view and Matplotlib for a pressure profile. Keep the numerical package usable directly from Julia without the Python UI.

The objective is a small, scientifically checkable vertical slice. PoreSpy, OpenPNM, IDAES, optimization, automatic differentiation, transient flow, real images, cloud hosting, authentication, and production job orchestration belong to later examples. The synthetic network needs no external dataset. This example validates an idealized model and the software connection; it does not establish predictive accuracy for a real porous material.

## 2. What has actually been checked

These observations apply **only to the hosted Linux workspace used to prepare this document**, not to the user's computer. The accessible workspace contained no project source files or Git repository.

| Tool or package | Observed status |
| --- | --- |
| Python | 3.12.14; executable available |
| uv | 0.13.5; executable available |
| Git | 2.51.1; executable available |
| Julia | Executable not found on PATH; installation elsewhere was not ruled out |
| NumPy | 2.3.5; import passed |
| SciPy | 1.17.0; import passed |
| VTK | 9.3.1; import passed; rendering not checked |
| pandas | 2.2.3; import passed |
| Matplotlib | 3.10.8; import passed |
| Panel, PyVista, JuliaCall, pytest | Not found in the active Python environment |
| PoreSpy, OpenPNM | Not found; not required for this example |
| Local Codex CLI | Executable not found on PATH in the hosted workspace |
| HTTPS package metadata | PyPI Panel metadata and Julia download metadata both returned HTTP 200; package installation was not tested |

An independent Python/SciPy sparse-network calculation was run during planning. The default case and six additional cases passed checks against analytical flow, pressure, conservation, and scaling. This was a reference check, **not a Julia solver, JuliaCall integration, UI, or installation test**. Its expected results are in section 5, so this handoff has no dependency on files from that workspace.

## 3. First action: local preflight and isolated setup

Inspect applicable AGENTS.md files before editing. Record the current directory, OS/architecture, relevant versions, whether this is a Git repository, and any existing changes. Use the intended empty project directory; if it contains unrelated work, put this example in a clearly named subdirectory without restructuring that work.

Create a standard-library-only `scripts/preflight.py` that can run before installing Python packages. It should report:

- Python executable and version; paths/versions of Julia, juliaup, uv, Git, and the local Codex CLI when discoverable.
- Relevant package versions and actual import success in the interpreter being tested. Report an absent package separately from a failed import.
- Whether the project environments and lockfiles exist.
- Optional bounded HTTPS checks for the package indexes/download endpoints actually needed. A successful metadata request does not prove that all package artifacts are reachable.

Do not dump the full environment, credentials, or unrelated user files. A missing optional tool is information, not an automatic blocker.

Prefer CPython 3.12 for the first example, a project-local `.venv` managed by uv, and one supported stable Julia release installed through the official platform-appropriate method. Record the exact Julia version actually used. The hosted package versions above are observations, not mandatory pins for the user's machine.

Use current official installation instructions. Reuse a suitable existing Julia installation. Install missing components in project/user scope where authorized and possible; avoid replacing system Python, privileged system changes, or modifying unrelated environments. If network or installation permissions block a required step, preserve completed work, state the exact blocker, and provide the minimum command needed to continue. Never claim a hybrid run passed by substituting a Python solver.

Direct Python runtime dependencies:

```text
numpy
panel
pyvista
matplotlib
juliacall
```

Use pytest as a development dependency. Let the resolver select compatible transitive packages; record them in `uv.lock`. SciPy is not needed at runtime: Julia owns the simulation. Julia's numerical core should initially need only `SparseArrays` and `LinearAlgebra`, plus `Test` for tests. The shared Julia project also needs `PythonCall` for the bridge.

### Use one explicit Julia environment

Put the Julia package/project in `julia/`. Both direct Julia execution and Python must use this same project. Before importing `juliacall`, configure:

- `PYTHON_JULIACALL_EXE`: the absolute path of the selected Julia executable.
- `PYTHON_JULIACALL_PROJECT`: the absolute path of the project's `julia/` directory.

The current JuliaCall documentation requires both variables when overriding its automatic environment management, and requires the Julia project's `PythonCall` version to match the Python `juliacall` version. Implement an idempotent `scripts/bootstrap_julia.py` that reads the installed Python package version, installs/checks the corresponding PythonCall version in this project, and verifies the selected environment. Consult the current documentation if these requirements change.

Avoid an accidental second Julia download or a hidden bridge environment. Keep the application bridge configuration in one module and set it before the first JuliaCall import. Do not rewrite the dependency manifests on every application launch. Preserve `uv.lock`, `julia/Project.toml`, and the generated `julia/Manifest.toml`; record the tested Python and Julia versions. Use a local Julia depot only if the machine's normal user depot is unavailable or the environment requires it, and configure it consistently.

## 4. Simulation specification

### Geometry

Generate a rectangular, nearest-neighbor cubic lattice. A node represents a pressure junction; every edge represents a cylindrical capillary. Junctions have no additional pressure loss or storage volume.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `nx, ny, nz` | `6, 4, 4` | Node counts along x, y, z |
| `spacing_m` | `1.0e-3 m` | Neighbor distance and capillary length, a |
| `radius_m` | `50.0e-6 m` | Uniform capillary radius, r |
| `viscosity_Pa_s` | `1.0e-3 Pa s` | Constant Newtonian viscosity, mu |
| `p_in_Pa` | `100 Pa` | Pressure on the x = 0 face |
| `p_out_Pa` | `0 Pa` | Pressure on the x = L face |

For zero-based logical coordinates `i, j, k`, use:

```text
x = i*a
y = (j + 0.5)*a
z = (k + 0.5)*a
```

Connect each node to its existing positive-x, positive-y, and positive-z neighbor exactly once. Use a documented deterministic node ordering. The default has 96 nodes and 224 edges. There are no periodic connections. The transverse boundaries are closed because no edges cross those boundaries.

Define sample length and representative cross-sectional area as:

```text
L = (nx - 1)*a
A = ny*nz*a^2
```

**The area convention matters.** Each transverse node represents a cell of area a^2, with a half-cell between the outer node centers and the boundary. Do not use `(ny - 1)*(nz - 1)*a^2`, which describes a different convention and fails for a single capillary chain. The declared transverse domain is `[0, ny*a]` by `[0, nz*a]`.

### Governing equations

Assume steady, incompressible, fully developed laminar flow in each capillary, with uniform viscosity, no gravity, and no junction/entrance losses. For an oriented throat from node i to node j:

```text
g_ij = pi*r_ij^4 / (8*mu*l_ij)       [m^3/(Pa s)]
q_ij = g_ij*(p_i - p_j)             [m^3/s]
sum_j g_ij*(p_i - p_j) = 0          at every interior node
```

Apply inlet and outlet Dirichlet pressures. Assemble the sparse conductance Laplacian and solve the reduced interior system:

```text
K_II*p_I = -K_IB*p_B
```

Handle `nx = 2` explicitly: every node is a boundary node, so there is no interior solve. Compute signed throat flows and nodal net outflows from the same orientation convention.

Define total inlet flow as the sum of nodal net outflows on the inlet face, and total outlet flow as minus that sum on the outlet face. This handles boundary-to-boundary transverse edges without double counting. For nonzero pressure drop:

```text
delta_p = p_in - p_out
k_eff = mu*L*Q_in / (A*delta_p)      [m^2]
```

At zero pressure drop, the solution is constant pressure with zero flow. Report the inferred permeability as unavailable rather than dividing by zero or presenting an analytically substituted value as a measured result.

Validate finite inputs, integer dimensions (`nx >= 2`, `ny >= 1`, `nz >= 1`), positive spacing/radius/viscosity, and `2*radius < spacing`. Reject invalid cases with actionable messages. Bound UI dimensions to `nx <= 20`, `ny <= 8`, `nz <= 8` for this prototype.

### Small process-scale connection

Use the computed permeability in a uniform Darcy column model, implemented as a separate Julia function:

```text
Q_column = k_eff*A_column*delta_p_column / (mu*L_column)
```

Default column length: 0.02 m; cross-sectional area: 1.0e-4 m^2; pressure drop: 1000 Pa; viscosity: the same as the network. This assumes the idealized repeating medium has the same effective permeability at column scale. Present it as a demonstration of upscaling, not a validated prediction for a real sample. If network permeability is unavailable, the column prediction is unavailable too.

## 5. Independent reference values and scientific tests

For the uniform lattice, each of the `ny*nz` x-directed paths is a series of `nx - 1` capillaries. Transverse pressures are equal, so transverse flow is zero analytically:

```text
p(x) = p_in - (p_in - p_out)*x/L
Q_exact = ny*nz*pi*r^4*(p_in - p_out) / (8*mu*L)
k_exact = pi*r^4 / (8*a^2)
```

Derive and document these expressions independently of the matrix assembly. The solver must solve the network equations rather than use these formulas to produce its answer.

| Default-case quantity | Expected value |
| --- | --- |
| Nodes / throats | 96 / 224 |
| Sample length | 0.005 m |
| Representative area | 1.6e-5 m^2 |
| Single-throat conductance | 2.4543692606170265e-12 m^3/(Pa s) |
| Pressure at successive x planes | 100, 80, 60, 40, 20, 0 Pa |
| Total inlet/outlet flow | 7.853981633974485e-10 m^3/s |
| Same flow in display units | 47.1238898038469 microlitres/min |
| Effective permeability | 2.4543692606170265e-12 m^2 |
| Default column flow | 1.227184630308513e-8 m^3/s |
| Same column flow in display units | 0.736310778185108 mL/min |

The planning reference obtained maximum pressure error about 5.0e-14 Pa and normalized conservation error about 1.2e-15. Those are observed reference values, not portable acceptance thresholds or a performance benchmark.

Required Julia numerical tests:

1. Default pressures and total flow match the analytical solution; use relative tolerance 1e-9 for nonzero analytical flow and permeability, and pressure error no larger than `1e-9*max(abs(delta_p), 1 Pa)`.
2. Conservation: both maximum absolute interior nodal imbalance and `abs(Q_in - Q_out)` are at most `1e-9*abs(Q_exact)` for nonzero-flow reference cases. Empty interior sets have zero maximum imbalance.
3. Doubling pressure drop doubles Q and leaves k unchanged. Adding the same constant to both boundary pressures leaves all flows unchanged. Reversing pressure drop reverses Q and leaves k positive.
4. Doubling radius multiplies both Q and k by 16. Doubling viscosity halves Q and leaves k unchanged.
5. Changing nx changes Q through L while leaving k unchanged at fixed spacing/radius. Include a `(2, 1, 1)` single-capillary case and an asymmetric `(9, 3, 2)` lattice.
6. Zero pressure drop produces constant pressure and absolute flows below 1e-20 m^3/s, with unavailable inferred k and column flow. No division by zero or NaN presented as a valid value.
7. Invalid and non-finite inputs fail cleanly. Confirm deterministic node/edge counts, valid connectivity, and near-zero transverse flow for a uniform case.

Required Python integration tests:

- The Python bridge obtains the same default answer as direct Julia execution.
- Returned array shapes, dtypes, units, edge orientation, and node IDs follow the documented contract.
- Exported summaries preserve the actual inputs and agree with the simulation values.
- UI-to-SI conversions are correct; in particular, micrometres to metres and m^3/s to the displayed flow units.

Do not add a second production solver in Python. These tests must call the actual Julia package; the analytical expressions provide the independent oracle.

## 6. Component boundary and proposed files

The exact internal organization can change if it improves simplicity, but preserve these responsibilities:

| Path | Responsibility |
| --- | --- |
| `README.md` | Setup, model assumptions, tested versions, launch/test commands |
| `pyproject.toml`, `uv.lock`, `.python-version` | Python environment and package configuration |
| `julia/Project.toml`, `julia/Manifest.toml` | Explicit Julia environment |
| `julia/src/Pore2ProcessDemo.jl` | Geometry, sparse solve, diagnostics, column calculation |
| `julia/test/runtests.jl` | Scientific tests |
| `src/pore2process_demo/bridge.py` | Runtime configuration, Julia loading, result conversion |
| `src/pore2process_demo/cli.py` | Run the default case and export outputs without a browser |
| `src/pore2process_demo/visualization.py` | Pressure plot and network visualization |
| `app.py` | Small Panel interface |
| `scripts/preflight.py` | Standard-library-only environment inventory |
| `scripts/bootstrap_julia.py` | Idempotent Julia project/bridge setup |
| `tests/` | Python bridge and export tests |
| `docs/environment_check.md` | Actual target-machine findings and setup changes |
| `docs/implementation_status.md` | Completed checks, results, limitations, next step |
| `results/default/` | Generated default-case outputs |

Ignore local environments, caches, compiled artifacts, and routine generated results in Git. Keep a small documented reference summary in the documentation. Use Git locally if available; remote hosting/pushing is not required.

Expose a coarse Julia API accepting scalar inputs, conceptually:

```text
run_case(nx, ny, nz, spacing_m, radius_m, viscosity_Pa_s, p_in_Pa, p_out_Pa)
predict_column_flow(permeability_m2, viscosity_Pa_s, length_m, area_m2, delta_p_Pa)
```

Return a structured result with coordinates `(N, 3)`, throat connectivity `(M, 2)`, node pressures `(N,)`, oriented throat flows `(M,)`, radius/length arrays, Q_in/Q_out, sample dimensions, effective permeability, conservation diagnostics, and timing information. Keep SI units at this boundary. Python converts display units and materializes ordinary NumPy arrays. Convert Julia's one-based node indices to Python/VTK zero-based indices **once in the bridge**, and test that conversion. No mutable simulation state should be shared between runs.

Use one persistent Julia runtime for the local application process. Begin with bounded, small synchronous solves on the calling server thread; avoid introducing JuliaCall background-thread/GIL/task-scheduler interactions in this first example. Keep controls/results per browser session. If observed execution time makes the interface unresponsive, document it and introduce a deliberately tested worker design; do not claim background execution or cancellation that has not been implemented.

## 7. Interface and exports

Build one page with:

- A compact explanation of the idealized capillary model and its units.
- Controls for lattice dimensions, spacing in mm, radius in micrometres, viscosity in mPa s, and inlet/outlet pressure in Pa. Defaults must reproduce section 5.
- An explicit **Run simulation** button. Editing inputs marks existing results as belonging to the previous run; it does not launch repeated solves automatically. Disable Run while it executes and restore it after success or failure.
- A pressure-colored 3D network with a labelled colorbar and readable axes. Use PyVista to build the geometry and Panel's documented VTK pane; trame is not required for this route. Use lines/glyphs or tubes with a clearly stated visual scale if exaggerated.
- A mean pressure versus x plot showing the numerical result and analytical line, with units.
- Total flow, permeability, conservation error, and default column flow, including unavailable states at zero pressure drop.
- A small secondary section for changing column length, area, and pressure drop using the last completed permeability. Record which network run supplied it. These controls should not rerun the network solver.
- A download/export action tied to the completed result and its input snapshot. Clearly show validation errors and preserve the last successful result without presenting it as the failed run.

Export at least:

| Output | Contents |
| --- | --- |
| `summary.json` | Input snapshot, explicit SI units, counts, flow, k, column inputs/result, diagnostics, runtime versions, timing, timestamp |
| `nodes.csv` | Zero-based node ID, coordinates in metres, pressure in Pa |
| `throats.csv` | Zero-based endpoint IDs, radius/length in metres, signed flow in m^3/s |
| `pressure_profile.png` | Labelled numerical/analytical comparison |
| `network.vtp` | Network geometry and pressure/flow fields for external inspection |

Represent unavailable values as JSON null, not nonstandard NaN/Infinity. File headers or metadata must make the flow orientation explicit. CLI exports can write `results/default/`; UI downloads should use a session-specific result snapshot and avoid a shared mutable output directory.

If graphics support is missing, complete the solver, bridge, CLI, exports, and 2D plot first. Diagnose the installed VTK/PyVista stack using its current platform guidance. A static fallback can keep the example usable, but interactive 3D remains an explicitly incomplete acceptance item until checked. Do not assume the hosted workspace's VTK version or rendering capability describes the local machine.

## 8. Execution sequence

Proceed in this order, resolving failures at the relevant step:

1. **Inventory:** inspect local instructions/files, create/run preflight, and document findings.
2. **Set up:** create the Python package/environment, install compatible dependencies, select Julia, configure the shared Julia project, and generate locks.
3. **Implement the Julia core:** build geometry and equations, run direct Julia tests, and compare to section 5 before UI work.
4. **Connect Python:** implement/configure the bridge, test array conversion and indexing, and run the complete default case through Python.
5. **Make it inspectable:** implement CLI, numerical exports, pressure plot, and the small Darcy column calculation.
6. **Build the interface:** add controls, Run, the network pane, result status, column controls, and downloads.
7. **Verify:** run the numerical/integration suite once the implementation is stable; exercise the UI through an available browser tool, including a changed input and exported result. If automated browser control is unavailable, document the exact manual checks still needed and perform HTTP/import checks without equating those with a visual pass.
8. **Handoff:** update README and implementation status with actual commands, versions, results, and remaining blockers. Do not end after scaffolding if required tools are available and the remaining work is authorized.

Target commands to support after implementation, from the project root (adapt shell quoting to the detected OS):

```bash
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run python scripts/preflight.py
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
uv run pytest
uv run python -m pore2process_demo.cli --output results/default
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

The direct `julia` command must resolve to the same executable selected for the bridge; otherwise document and use its explicit path. A first-time installation may require initial `uv add`/lock creation and a Julia installer before these commands are available. These commands describe the desired completed project, not files already delivered with this plan. Record cold initialization separately from warmed solver time; no speedup claim is part of this example.

## 9. Definition of done

- [ ] Target-machine inventory is recorded, and missing tools are installed or precisely identified as blockers.
- [ ] Python and Julia environments are isolated/explicit, compatible, and reproducible through the saved manifests/locks.
- [ ] Direct Julia numerical tests pass against the independent analytical expressions.
- [ ] Python calls the real Julia package and receives correctly converted results.
- [ ] The default case reproduces section 5 within the stated tolerances.
- [ ] A changed pressure, radius, or viscosity produces the expected behavior.
- [ ] The effective permeability drives the demonstrated Darcy column calculation.
- [ ] CLI output is reproducible and includes the required exports and units.
- [ ] The local Panel app runs, displays the 3D network and pressure profile, handles invalid/zero-flow cases, and exports the completed result.
- [ ] Documentation distinguishes numerical tests, browser checks, and anything not verified.
- [ ] The user can reproduce setup, run, and tests using the actual documented commands.

The next example, after this passes, should introduce one realistic complication such as heterogeneous throat radii or imported pore-network geometry. Choose that extension from the research question rather than expanding this first example prematurely.

## 10. Official references

Consult current primary documentation for package APIs and installation details; do not treat remembered snippets as authoritative when they fail.

- Julia installation: https://julialang.org/install/
- Julia package environments: https://pkgdocs.julialang.org/v1/environments/
- Julia sparse arrays: https://docs.julialang.org/en/v1/stdlib/SparseArrays/
- JuliaCall setup, environment overrides, and threading: https://juliapy.github.io/PythonCall.jl/stable/juliacall/
- uv projects and locks: https://docs.astral.sh/uv/guides/projects/
- Panel: https://panel.holoviz.org/
- Panel VTK pane: https://panel.holoviz.org/reference/panes/VTK.html
- PyVista: https://docs.pyvista.org/
- PyVista/trame documentation, if that alternative is later needed: https://docs.pyvista.org/api/plotting/trame.html
- OpenPNM/PoreSpy integration, for a later example: https://openpnm.org/examples/applications/network_extraction.html
- IDAES overview, for a later process integration: https://idaes.org/overview/
