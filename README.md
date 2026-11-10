<p align="center">
  <img src="assets/branding/pore2process-logo.png" alt="Pore2Process — packed porous column and wordmark" width="820">
</p>

# Pore2Process

**From porous-material structure to process-scale performance.**

Pore2Process is a research software project for connecting what happens inside a
porous material to what matters in a chemical process: pressure requirements,
pumping energy, throughput and, in future applications, conversion and thermal
performance.

The aim is a reproducible chain of calculations in which a change in pore
geometry produces an explainable change in a process result, with the
assumptions and intermediate evidence available for inspection.

> **Project status — 10 November 2026.** The filtration hydraulic example is
> complete and closed for now, with source, tests, dependency locks and curated
> references preserved as a regression benchmark. **Isothermal adsorption** is
> the next case in development. The verified native setup targets Apple Silicon
> macOS.

## The connected workflow

```mermaid
flowchart LR
    G["Synthetic or literature-informed geometry"] --> O["OpenLB<br/>Pore-scale simulation"]
    O --> J["Julia<br/>Effective properties and scale translation"]
    J --> D["DWSIM<br/>Process simulation"]
    D --> R["Process metrics<br/>Energy, throughput and performance"]
    P["Python<br/>Coordination, data and interface"] -.-> O
    P -.-> J
    P -.-> D
```

| Component | Defined role |
| --- | --- |
| **OpenLB** | Detailed pore-scale physics using existing OpenLB facilities and the structure of its official examples. |
| **Julia** | Custom numerical models, analysis, effective-property extraction and translation between scales. |
| **Python** | Workflow coordination, data handling, integration, visualization and a local Panel interface. |
| **DWSIM** | Actual process calculations in saved flowsheets, including equipment duties, recycle and balances. |

Models remain usable independently of the interface.

## Examples and scientific reports

Each example has a `REPORT.md` describing its purpose, assumptions, theory and
equations, recorded results, discussion and conclusions. Setup and launch
commands remain in the corresponding `README.md`.

| Area | Starting points |
| --- | --- |
| **Hydraulics** | [Capillary network](examples/hydraulics/capillary_network/REPORT.md), [3D capillary verification](examples/hydraulics/capillary_benchmark/REPORT.md), [measured permeability to pump duty](examples/hydraulics/darcy_pump/REPORT.md) |
| **Filtration** | [Study stages](examples/filtration/README.md), [fixed-opening support design](examples/filtration/graded_support/REPORT.md), [final operating map and conclusions](examples/filtration/operating_map/REPORT.md) |
| **Adsorption** | [Synthetic bed](examples/adsorption/REPORT.md), [literature LDF comparison](examples/adsorption/alhamed/REPORT.md), [film/particle transport](examples/adsorption/film_particle/REPORT.md), [column contact and inventory](examples/adsorption/column_contact/REPORT.md) |

The [complete example catalogue](examples/README.md) links all scientific
reports and run guides. Reports distinguish numerical verification,
literature comparisons and conditional scenarios; unfinished and failed
comparisons remain explicit.

## How results stay traceable

- **Explicit contracts:** SI units, coordinate systems, node indexing, array
  shapes, flow orientation and file formats are recorded at software boundaries.
- **Independent checks:** analytical solutions, conservation, numerical
  convergence and appropriate published comparisons support each stage.
- **Saved evidence:** inputs, assumptions, solver settings, diagnostics, software
  versions and outputs travel with each completed run.
- **Clear provenance:** extracted measurements are distinguished from illustrative
  scenarios; unreported quantities remain unavailable or are explored as labelled
  assumptions.
- **Reproducible environments:** the implementation uses an isolated Python
  environment, a shared explicit Julia project and saved dependency locks.

Published articles and their available supplements provide experimental
evidence.

## Scientific references

The membrane-structure, hydraulic and crossflow studies draw on the following
articles. Source locations and extraction records are in the
[literature workspace](literature/README.md) and [source ledger](literature/sources.json).

1. **Biesheuvel, P. M., and Verweij, H. (1999).** [Design of ceramic membrane supports: permeability, tensile strength and stress](https://doi.org/10.1016/S0376-7388(98)00335-4). *Journal of Membrane Science*, 156, 141–152.
2. **Zhu, W., et al. (2021).** [Design and optimization of ceramic membrane structure: From the perspective of flux matching between support and membrane](https://doi.org/10.1016/j.ceramint.2021.01.088). *Ceramics International*, 47, 12357–12365.
3. **Mo, J., Li, X., and Yang, Z. (2022).** [Dissecting the structure-property relationship of ceramic membrane with asymmetric multilayer structures for maximizing permselectivity](https://doi.org/10.1016/j.watres.2022.118658). *Water Research*, 220, 118658.
4. **Jan, A., et al. (2025).** [Silicon carbide coated alumina tight-ultrafiltration membrane prepared by low-pressure chemical vapor deposition for sulphate ion retention](https://doi.org/10.1016/j.desal.2025.119085). *Desalination*, 613, 119085.
5. **Chiu, T. Y., and James, A. E. (2005).** [Critical flux determination of non-circular multi-channel ceramic membranes using TiO2 suspensions](https://doi.org/10.1016/j.memsci.2005.01.037). *Journal of Membrane Science*, 254, 295–301.
6. **Chiu, T. Y., Lara Dominguez, M. V., and James, A. E. (2006).** [Critical flux and rejection behaviour of non-circular-channelled membranes: Influence of some operating conditions](https://doi.org/10.1016/j.seppur.2005.11.028). *Separation and Purification Technology*, 50, 212–219.
7. **Garcia Garcia, F. J., and Chiu, T. Y. (2008).** [Economic aspects of critical flux operability in star shaped microfiltration membranes: Influence of some operating conditions](https://doi.org/10.1016/j.memsci.2008.08.033). *Journal of Membrane Science*, 325, 641–646.
8. **Frederic, E., et al. (2018).** [Modelling of fluid flow distribution in multichannel ceramic membrane: Application to the filtration of produced water](https://doi.org/10.1016/j.memsci.2018.09.021). *Journal of Membrane Science*, 567, 290–302.

## Where this can go

The next case adds **isothermal adsorption and mass transfer**. Its research
question is how particle size changes usable treatment capacity and pumping
energy at fixed adsorbent inventory and required outlet quality. That work is
in progress; its current source and comparison reports are in
[examples/adsorption/](examples/adsorption/README.md).

The filtration example remains available for regression and comparison.
Reacting and non-isothermal chemical and bioprocess applications remain
longer-term directions.

See the [project roadmap](ROADMAP.md) for milestones and acceptance criteria.

## Project layout

Runnable studies are grouped in [examples/](examples/README.md):

- **[Hydraulics](examples/hydraulics/capillary_network/README.md):** network reference, native capillary benchmark, Darcy/pump coupling and tool checks.
- **[Filtration](examples/filtration/README.md):** the completed support/channel/crossflow operating-map study.
- **[Adsorption](examples/adsorption/README.md):** bed hydraulics, uptake models and literature comparisons in development.

Shared Python code stays in `src/pore2process_demo/`; Julia models and their
explicit environment stay in `julia/`. The article workspace, existing runs and
local tool installations remain in `literature/`, `results/` and `.local/`.
See [the layout and migration guide](docs/project_layout.md).

## Start with the capillary-network example

Install **uv** and **Julia** using their official installation instructions,
then run from a terminal. The initial Python–Julia example does not require
OpenLB or DWSIM.

```sh
git clone https://github.com/axeld92/Pore2Process.git
cd Pore2Process
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run python -m pore2process_demo.cli --output results/default
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

Open `http://127.0.0.1:5006/app` and press **Run simulation** in **Capillary
network**. Bootstrap accepts `--julia /absolute/path/to/julia` if Julia is not on
PATH. Python 3.12 and Julia 1.11.6 were tested with matching
JuliaCall/PythonCall 0.9.37; the exact dependencies are saved in `uv.lock` and
`julia/Manifest.toml`.

```sh
uv run python scripts/preflight.py
julia --startup-file=no --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
uv run pytest -q
```

Use the same Julia executable selected by bootstrap for direct Julia tests.
The VTK tests need normal host graphics access. The recorded full numerical
and interface checks are in the [acceptance report](docs/implementation_status.md)
and [browser record](docs/browser_checks.md).

## Reproduce the connected showcase

The native setup has been checked on **macOS ARM64**, with **OpenLB 1.9r0**,
**DWSIM 10.2.15**, a private **.NET 10** runtime and **Homebrew GCC 15**. Its
builders currently use that platform's compiler and library paths.

Follow the [clean-checkout reproduction guide](docs/reproduction.md) to download
the pinned tools, run the prerequisite studies and generate the filtration
showcase. The interface uses completed native evidence and launches actual
Julia/DWSIM process calculations; it does not rerun OpenLB whenever a control
changes. Native study generation can take substantial time.

| Documentation | Contents |
| --- | --- |
| [Setup and examples](docs/setup_and_examples.md) | Detailed Python–Julia setup, model derivation, contracts and original example commands |
| [Clean-checkout reproduction](docs/reproduction.md) | Native prerequisites and ordered connected-study commands |
| [Filtration showcase](examples/filtration/star_channel/REPORT.md) | Geometry, equations, results and biological-feed assumptions |
| [Closed filtration operating map](examples/filtration/operating_map/REPORT.md) | Final hydraulic range, material/operation comparison and circulation-energy result |
| [Literature workspace](literature/README.md) | Curated references, extractions and provenance |
| [Integration roadmap](docs/integration_roadmap.md) | Detailed milestone and research status |
| [Publication record](docs/code_publication.md) | Published snapshot and verification scope |
