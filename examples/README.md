# Examples and scientific reports

Each `REPORT.md` explains the question or purpose, model assumptions, theory
and equations, recorded results, discussion and conclusions. The reports use
retained evidence; preparing them does not rerun the studies. A software or
numerical pass is distinguished from a successful literature prediction.

Use each `README.md` for setup and launch/check commands, from the project root.
The reports can be read individually; links explain how later examples use
specific earlier measurements or assumptions.

| Area | Scientific report | Purpose / finding | Run guide |
| --- | --- | --- | --- |
| Hydraulics | [Capillary network](hydraulics/capillary_network/REPORT.md) | Julia flow reference and Darcy upscaling | [Commands](hydraulics/capillary_network/README.md) |
| Hydraulics | [3D capillary](hydraulics/capillary_benchmark/REPORT.md) | Native measured-permeability verification; current TRT evidence | [Commands](hydraulics/capillary_benchmark/README.md) |
| Hydraulics | [Darcy and pump](hydraulics/darcy_pump/REPORT.md) | Pore-derived properties drive actual DWSIM work | [Commands](hydraulics/darcy_pump/README.md) |
| Hydraulics | [Tool checks](hydraulics/tool_checks/REPORT.md) | Independent native channel and pump compatibility | [Commands](hydraulics/tool_checks/README.md) |
| Filtration | [Layers in series](filtration/layered_support/REPORT.md) | Support benefit and fixed-layer dominance | [Commands](filtration/layered_support/README.md) |
| Filtration | [Literature calibration](filtration/calibration/REPORT.md) | Failed support closure; bounded matched-layer reference | [Commands](filtration/calibration/README.md) |
| Filtration | [Fixed-opening support](filtration/graded_support/REPORT.md) | Controlled interior geometry and common scale translation | [Commands](filtration/graded_support/README.md) |
| Filtration | [Once-through circuit](filtration/once_through/REPORT.md) | Recovery, feed duty and prescribed channel loss | [Commands](filtration/once_through/README.md) |
| Filtration | [Crossflow recycle](filtration/crossflow_recycle/REPORT.md) | Separate circulation energy and actual recycle balances | [Commands](filtration/crossflow_recycle/README.md) |
| Filtration | [Literature capacity](filtration/literature_benchmark/REPORT.md) | Critical-flux policies and unresolved source-flow basis | [Commands](filtration/literature_benchmark/README.md) |
| Filtration | [Explicit low-Re star](filtration/star_channel/REPORT.md) | Distributed hydraulics and illustrative COD balance | [Commands](filtration/star_channel/README.md) |
| Filtration | [Final operating map](filtration/operating_map/REPORT.md) | Higher-Re verification and conditional process selections | [Commands](filtration/operating_map/README.md) |
| Adsorption | [Synthetic bed](adsorption/REPORT.md) | Connected uptake/hydraulics/pump scenario | [Commands](adsorption/README.md) |
| Adsorption | [Date-stone LDF](adsorption/alhamed/REPORT.md) | Failed calibration and withheld-size column prediction | [Commands](adsorption/alhamed/README.md) |
| Adsorption | [Film and particle](adsorption/film_particle/REPORT.md) | Withheld batch loading passes; column transfer fails | [Commands](adsorption/film_particle/README.md) |
| Adsorption | [Column contact](adsorption/column_contact/REPORT.md) | Common mixing and inventory audit; column comparison fails | [Commands](adsorption/column_contact/README.md) |

Filtration hydraulics is closed for now; adsorption remains research in
progress. Earlier BGK numbers and later TRT results retain separate identities.
The reports give their own experimental/scenario and operating-range limits.

Shared numerical models remain in `julia/`; shared Python runtime and interface
code remains in `src/pore2process_demo/`. Existing run folders, literature
extractions and tool installations retain their locations. See the
[layout guide](../docs/project_layout.md) and
[ordered native/process reproduction](../docs/reproduction.md).
