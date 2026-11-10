# Examples

Run commands from the project root after `uv sync --locked` and the shared Julia
bootstrap. All native simulator cases keep their code and build outputs outside
the original OpenLB distribution. Installed tools stay in `.local/` and completed
runs stay under their existing `results/` paths.

| Area | Example | Purpose / status |
| --- | --- | --- |
| Hydraulics | [Capillary network](hydraulics/capillary_network/README.md) | Original Python–Julia numerical/UI reference |
| Hydraulics | [Native capillary](hydraulics/capillary_benchmark/README.md) | Native pressure-flow and refinement benchmark |
| Hydraulics | [Darcy pump](hydraulics/darcy_pump/README.md) | Measured permeability drives actual DWSIM pump work |
| Hydraulics | [Tool checks](hydraulics/tool_checks/README.md) | Independent OpenLB/DWSIM installation probes |
| Filtration | [Overview](filtration/README.md) | Closed hydraulic study and ordered reproduction |
| Adsorption | [Overview](adsorption/README.md) | Synthetic bed and literature comparisons; research work in progress |

Shared numerical models remain in `julia/`. Runtime setup, JSON/hash helpers,
external builders and the local interface remain in `src/pore2process_demo/`.
`app.py` remains the local Panel entry point. See [the layout guide](../docs/project_layout.md).
