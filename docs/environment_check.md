# Target-machine inventory — 7 November 2026

## Crossflow recycle continuation — 9 November 2026

Root instructions and Git status were rechecked. The only existing untracked
directory was unrelated `assets/`, which remains untouched. The actual target
still reports macOS 26.6.2 ARM64, project CPython 3.12.9, Julia 1.11.6, private
.NET SDK 10.0.403 / runtime 10.0.14 and DWSIM FluentAPI 10.2.15.0.
No environment, package, lock or original simulator-source changes were needed.

The installed recycle/mixer/cooler APIs and bundled official recycle example
were inspected before implementation. Reflection inventories, build outputs
and exploratory pilots are isolated under `.local/crossflow_recycle/`. A stock
pure-water tear-stream flash setting resolves the observed temperature-residual
stall while retaining the original convergence tolerances; details and primary
references are in [the result report](crossflow_recycle.md). The standard Julia
user-cache and original VTK graphics sandbox requirements are unchanged.

## Explicit DWSIM filtration — 8–9 November 2026

Root `AGENTS.md`, workspace and Git status were inspected before editing on
`codex/first-simulation`. The only pre-existing untracked directory was the
unrelated `assets/`, which was preserved. The target remains macOS 26.6.2
(Darwin ARM64), CPython 3.12.9 and the same Julia 1.11.6 shared project.
The existing private SDK reports 10.0.403; the installed DWSIM FluentAPI reports
10.2.15.0. Its supported builders were inspected by reflection before selecting
stock operations; [the interface decision](dwsim_filtration_interfaces.md)
records details and current official references.

No packages, runtimes or locks were changed. New adapter builds are isolated
under `.local/process_filtration/`; the original adapter and vendor sources
are unchanged. Julia's Python-call precompilation again required normal user
cache access; it passed after that sandbox restriction was lifted. The full
Python/VTK suite uses normal macOS graphics access. Those automation sandbox
requirements are distinct from an application installation blocker.

## Filtration continuation — 8 November 2026

The workspace was inspected again on `codex/first-simulation` with a clean
working tree before this continuation. Root `AGENTS.md` now supplies the OpenLB
policy and was read; the original inventory below predates those instructions.
The target remains macOS 26.6.2 ARM64, project Python 3.12.9, uv 0.13.7 and the
same Julia 1.11.6. Fresh preflight output is retained at
`results/filtration/preflight.json`. Existing installed dependencies and private
OpenLB/DWSIM runtimes were reused; no software installation was needed.

The filtration case compiles in `.local/filtration/build/` using the shared
external builder and existing application header. It does not edit the vendor
tree. The pinned-archive comparison checks 2,070 files, preserving the existing
documented exception for a vendor config file changed during an older probe.

The first sandboxed launch reached DWSIM but was prevented from creating a
Julia precompile pidfile in the normal user depot. Running with normal user
cache access resolved that restriction. An unsourced uv invocation likewise
could not create the default user cache; sourcing `scripts/env.sh` uses the
project cache. These were execution-environment restrictions, not dependency
or solver failures. The original PyVista regression tests still require normal
macOS graphics access. No new connected interface is introduced by this study.

## Original setup inventory

The target was `/Users/dula/Documents/ChatGPT/Pore2Process` on macOS 26.6.2,
Darwin arm64. It was a Git repository on unborn `main` with no commits or
existing application code. The supplied handoff was read in full and preserved.
No applicable AGENTS.md was found in the project or its ancestor directories.
The hosted Linux observations were not used as target-machine facts.

Initial shell Python resolved through pyenv to 3.11.9. After `.python-version`
requested the then-unavailable 3.12, an unqualified preflight invocation ran
Apple's command-line-tools Python 3.9.6; `docs/preflight_initial.json` explicitly
records that executable. Its available NumPy/Matplotlib are observations of that
interpreter only. uv subsequently discovered existing Homebrew Python 3.12.9
at `/opt/homebrew/opt/python@3.12/bin/python3.12` and created this project's `.venv`.
No Python runtime had to be downloaded.

The existing Juliaup launcher resolved to Julia 1.11.6 at
`/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia`.
Juliaup is 1.24.0. Git is 2.49.0 and the discoverable local Codex CLI is 0.160.1.
uv was absent from PATH and checked common locations. It was installed from
PyPI into ignored project-local `.tools/` (version 0.13.7), using the existing
pyenv Python 3.11.9 with `pip install --target .tools uv`.

Python runtime dependencies and pytest were installed into `.venv` and resolved
in `uv.lock`. Actual imports passed for NumPy 2.5.4, Panel 1.10.0, PyVista
0.49.2, VTK 9.7.2, Matplotlib 3.11.2, JuliaCall 0.9.37 and pytest 9.2.0.
SciPy is not a runtime dependency. Panel resolved Bokeh 3.10.2. The full
transitive set is in the lockfile, rather than duplicated here.

Bootstrap added PythonCall 0.9.37 to the explicit `julia/` project and saved
its exact compatibility requirement and generated manifest. Both direct Julia
and the Python bridge use this project and the same executable. The normal
`~/.julia` depot was available with installation permission and was reused.
PythonCall's CondaPkg dependencies fetched their bundled pixi/micromamba
artifacts during precompilation, but `JULIA_CONDAPKG_BACKEND=Null` and the
explicit Python executable prevent provisioning another Python environment.
JuliaCall's automatic Julia/environment management is disabled by setting both
documented overrides. Bootstrap was rerun successfully and manifest hashes
were unchanged; subsequent launches do not invoke Pkg.

The terminal sandbox initially denied package-network access, Juliaup/cache
writes, and macOS graphics initialization. Package installation and full UI
tests succeeded with approved normal host access. No system Python replacement,
privileged installation, or shell-profile edits were performed. Matplotlib and
font caches are redirected into `.local/` by the environment script/application.
The raw preflight intentionally reports its interpreter's import behavior;
without the environment script it can report a temporary Matplotlib cache warning.

The first browser run exposed Bokeh 3.10 moving synchronous callbacks onto worker
threads. A second Julia run deadlocked in Julia GC while the main Python thread
waited for the GIL. A short native process sample established those call stacks.
Event callbacks now use async wrappers executed on the main event loop, with
synchronous Julia calls inside them, and the bridge rejects calls from other
threads. The final acceptance report distinguishes the repeated browser checks
from direct method/unit tests.

`docs/preflight_verified.json` records the installed project-interpreter imports.
Network artifacts were actually installed; HTTPS metadata alone was not treated
as installation proof. The app's HTTP response was independently checked:
status 200 with an HTML document. Native graphics and browser rendering were
checked separately; HTTP success was not used as a visual pass.


## 8 November: external OpenLB builds and connected demonstration

The machine inventory was repeated with the project interpreter. macOS26.6.2
ARM64, Python3.12.9, Julia1.11.6, JuliaCall/PythonCall0.9.37, uv0.13.7,
Homebrew GCC15.2.0, tinyxml2 11.0.0 and zlib1.2.12 remain the tested configuration.
Private DWSIM10.2.15 and .NET SDK10.0.403/runtime10.0.14 are reused. Locked Python
synchronization passes; no dependency or system-installation changes were needed.

Fresh inventories are in `results/preflight_coupling.txt` and
`results/compatibility/final_environment.txt`. Build configuration, objects,
core archive, executable and copied stock example now live outside the original
OpenLB distribution. The preexisting vendor config from the older compatibility
probe is ignored. All 2,070 checked original distribution files match the pinned
archive. The project code policy and exact setup/run commands are documented
in `docs/openlb_code_policy.md`, `examples/hydraulics/capillary_benchmark/README.md` and `examples/hydraulics/darcy_pump/README.md`.

Julia's normal user compilation cache required approved host access to update
a compiled project module; the original VTK regression suite required normal
macOS graphics access. Those checks passed. No genuine installation, build,
solver or coupling blocker remains for the documented local demonstration.

## Graded-support setup rechecked — 8 November 2026

The local preflight was rerun before this extension. The machine remains
macOS 26.6.2 / ARM64, CPython 3.12.9, Julia 1.11.6 and project-local uv 0.13.7.
All required Python imports pass and the existing locks/shared Julia project
are present. The actual DWSIM anchor reports version 10.2.15.0; the private
adapter build checks SDK 10.0.403. GCC 15.2 and the existing OpenLB 1.9r0 source
are reused without installation or dependency changes.

The new application code is under `examples/filtration/graded_support/openlb/`; binaries and
objects stay in `.local/graded_support/build/`. Vendor verification again
checks 2,070 original files against the pinned archive. The original numerical
cases, bridge and DWSIM adapter are preserved. The shared native launch helper
adds only an optional execution timeout, retaining its 600-second default for
earlier examples. This study requests 1,800 seconds per native solve after its
finest baseline took 535 seconds; numerical settings/screens remain fixed.
Two external native processes may run concurrently; Julia analysis runs only
on the main Python thread. See [the study commands](../examples/filtration/graded_support/README.md)
and [scientific scope](graded_support_study.md).

## TRT and explicit-channel extension — 9 November 2026

`results/star_channel/preflight_host.json` records a fresh project-interpreter
inventory: macOS 26.6.2 ARM64, Python 3.12.9, Julia 1.11.6, Juliaup 1.24.0,
uv 0.13.7, Git 2.49.0 and local Codex CLI 0.160.1. All required package imports
pass. NumPy 2.5.4, Panel 1.10.0, PyVista 0.49.2, VTK 9.7.2, Matplotlib 3.11.2,
JuliaCall/PythonCall 0.9.37 and pytest 9.2.0 remain the tested versions.
The existing GCC 15.2.0, OpenLB 1.9r0, DWSIM 10.2.15.0 and .NET SDK 10.0.403
are reused. No dependency installation, system change or lockfile update was
needed. Julia still uses the same explicit project and normal user depot.

The user requires at least TRT for porous media. Existing pressure-driven
application cases now use stock TRT with magic parameter 3/16. The new star
channel uses stock TRT collision with stock exact-difference forcing, composed
through OpenLB's existing `Tuple`. Project source/configuration is in
`examples/filtration/star_channel/openlb/`; its external build is `.local/star_channel/build/`.
All 2,070 original checked OpenLB files still match their pinned archive.

The eight-case capillary TRT regression and its independent SI/export replay
pass under `results/hydraulic_trt/`. Its generic launcher now permits 1,800
seconds per native solve. The separate graded-support launcher permits 3,600
seconds after its finest 6 µm case exceeded 1,800 seconds while several studies
were running. Completed native cases are reused only with matching binary,
input and artifact hashes; analysis and process solves are repeated. These
are execution limits, not steady-state or scientific tolerance changes.
Earlier BGK runs retain their original identities and values.

Normal approved host execution is needed for Julia's compilation cache,
native macOS graphics, and the localhost Panel server. Sandbox permission
limitations are distinguished from numerical failures. Native timing includes
overlapping runs and is not used for a performance comparison. For example,
the accepted 64-cell low-gradient star case recorded 43.49 s and its Re≈20
operating case 121.39 s; the separate cold Julia initialization and each
solver's timing remain in the generated summaries.

The complete TRT support study passes all 24 native cases; the finest 6 µm
case recorded 2,011 s after resuming with hash-verified earlier measurements.
The new process matrix passes 27 design/sensitivity cases plus seed validation.
It opts into DWSIM's stock forced-liquid route, checks saturation/freezing
margins and saved/reloaded phase settings, and tightens the recycle mass
tolerance to 1e-14 kg/s for low-duty accuracy. Older cases retain their default
phase calculation and original tolerances. No simulator source is patched.
The final regression suite passes 61 Python tests and 878 Julia assertions.
The live server uses the documented `uv run panel serve app.py --address
127.0.0.1 --port 5006` command, with browser evidence recorded separately.

## Higher-Re operating-map continuation — 9–10 November 2026

`results/operating_map/preflight.json` records the actual target-machine check.
The existing macOS ARM64 Python/Julia/OpenLB/DWSIM runtimes and dependency locks
are reused; no package or system installation was performed. All new native
application code is under `examples/filtration/operating_map/openlb/`, with external build products
in `.local/operating_map/build_trt27/`. Stock D3Q27/TRT replaces the rejected
higher-Re D3Q19 pilots for this new study; completed earlier examples retain
their own evidence and controls. Original OpenLB source remains untouched.

Native cases may use four independent single-core processes. The actual DWSIM
matrix uses two independent child processes. Julia calls stay on the main Python
thread. The initial native pass permitted 3,600 seconds per case; the longer-domain
case was restarted from its latest native field with a 7,200-second allowance.
Scientific stationarity, conservation, Mach, refinement and physical-time screens are
distinct from those execution ceilings. Normal user Julia-cache/graphics access is still
needed for the existing local runtime tests; no new installation blocker was
found. Unaccepted exploratory attempts remain under `results/operating_map/`.

All 16 native cases subsequently passed 144 acceptance checks. The restarted
domain case converged after 833,536 steps and 3,545.46 seconds in that attempt;
the earlier attempt and initial-field provenance remain saved. These times
include the actual execution conditions and are not a performance benchmark.

The first process launch stopped before a solve when the installed .NET SDK
version probe exceeded its existing 15-second ceiling. A normal-host retry
reported the unchanged SDK 10.0.403 in 0.35 seconds. The same process run was
resumed, with no SDK, simulator, scientific tolerance or source-code change.
The subsequent first process requests exposed our adapter's older 5,000-
iteration input guard. It now permits up to 20,000, as requested by the new
high-circulation protocol; older input defaults and the stock solver remain
unchanged. The failed attempt is preserved and the corrected process matrix
uses a fresh source/input snapshot. No result from the rejected attempt is
promoted as accepted evidence.

The direct numerical suite, including the bounded source relation and an
independent circular-profile momentum screen, passes 902 Julia assertions.
All 90 process design cases plus the independent seed case subsequently pass
actual DWSIM solve/save/reload checks; a one-iteration recycle is rejected.
Full independent replay passes, and the normal-host Python suite passes
74 tests. Native/process acceptance results are recorded in
[operating_window.md](operating_window.md). No dependency lock changed.

## Isothermal adsorption starting example — 10 November 2026

`results/adsorption/preflight.json` records a fresh target-machine inventory.
All existing Python/graphics/Julia imports pass. macOS ARM64, Python 3.12.9,
Julia 1.11.6, JuliaCall/PythonCall 0.9.37, GCC 15.2, OpenLB 1.9r0,
DWSIM 10.2.15 and the project-local .NET SDK 10.0.403 are reused. No Python
or Julia dependency installation, manifest change or privileged setup was
required. New C# restore/build products remain in `.local/adsorption/`.

The new applications are in `examples/adsorption/openlb/` and `examples/adsorption/openlb_bed/`;
the existing external builder puts all core/application build products in
separate `.local/adsorption/` folders. The first scalar compilation recorded
515.65 seconds; subsequent cached checks take much less. Native timings are
execution records, not performance comparisons. All 2,070 checked original
OpenLB files still match the pinned archive, with the pre-existing vendor
configuration exception retained.

The initial FCC-bed scaling failed its density-departure screen. It is
preserved in `results/adsorption/bed_initial_compressibility_rejected/`.
Smaller lattice time steps, with stock TRT and magic 3/16, pass the unchanged
1e-4 density limit. Optional initial velocity fields from that attempt are
only starting estimates; new convergence and a separate cold-start check
are required. Eight revised bed cases and seven scalar-TRT cases pass.

The new standalone Julia module stays in `julia/adsorption/`; loading it does
not edit the frozen filtration package or shared environment. Its 46 assertions
pass, including a fast-uptake zero-tail regression after improving the
conservative implicit line search. The existing package still passes all
902 assertions. Full normal-host Python regression passes 86 tests. No new
Julia worker thread or background-runtime interaction is introduced.

The low-head DWSIM diagnostic and stock forced-liquid remedy are documented in
[the adsorption report](adsorption_example.md). The separate adapter preserves
the earlier pump/recycle adapters and checks phase margins and reload work.
Four actual solves/reloads and independent replay pass. All native fields,
input/source snapshots and actual simulator archives are saved in a unique
accepted process run. Existing simulator distributions remain untouched.

The publicly available author PDF was saved in the literature inbox and its
tables visually checked. Bundled Poppler initially retried an unavailable
font cache; a project-local fontconfig file/cache resolved rendering without
system installation. This reading-tool issue is separate from simulation.
The scientific figure was visually inspected. There is no new adsorption
Panel interface or browser-pass claim.
