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
| Adsorption | [Native packed-bed tracer](adsorption/packed_bed_tracer/REPORT.md) | One pore-informed contact coefficient predicts a longer segment while preserving flow and liquid inventory; bounded synthetic regime | [Commands](adsorption/packed_bed_tracer/README.md) |
| Adsorption | [Column transport and pump duties](adsorption/column_transport/REPORT.md) | Passing outlet controls; failed application preservation blocks adsorption; independent fixed-flow hydraulic process branch | [Commands](adsorption/column_transport/README.md) |
| Adsorption | [Engineering transport accuracy](adsorption/transport_accuracy/REPORT.md) | Uniform drift passes the new engineering screen; fronts and a periodic porous profile reveal larger transport errors | [Commands](adsorption/transport_accuracy/README.md) |
| Adsorption | [Transport resolution and root cause](adsorption/transport_root_cause/REPORT.md) | Simple controls converge; matched N96 reduces the largest bed error but still fails; checked stock wall alternative rejected | [Commands](adsorption/transport_root_cause/README.md) |
| Adsorption | [Bulk stretching and one curved wall](adsorption/transport_verification/REPORT.md) | Independent references show late wall-free stretching error and separate distortion of uniform tracer during circular-wall rotation | [Commands](adsorption/transport_verification/README.md) |
| Adsorption | [Stationary curved wall](adsorption/stationary_curved_wall/REPORT.md) | Zero wall velocity reduces uniform drift by 86%; independent radial comparison still rejects the varying tracer | [Commands](adsorption/stationary_curved_wall/README.md) |
| Adsorption | [Path-specific wall uptake](adsorption/wall_uptake/REPORT.md) | Perfect-sink film estimates from saved TRT flow, with failed two-period transfer and unfitted Figure 7 checks | [Commands](adsorption/wall_uptake/README.md) |
| Adsorption | [Wall exposure among survivors](adsorption/survivor_exposure/REPORT.md) | Frozen uptake model passes the two-cell holdout; both experimental breakthrough comparisons still fail | [Commands](adsorption/survivor_exposure/README.md) |
| Adsorption | [Fast-path bypass across eight cells](adsorption/repeated_geometry/REPORT.md) | Frozen tracer model passes; within-class reset slightly increases predicted wall-free passage at fixed exchange | [Commands](adsorption/repeated_geometry/README.md) |
| Adsorption | [Ordered and disordered packing](adsorption/disordered_packing/REPORT.md) | Stock TRT comparison at fixed carbon and Péclet numbers; disorder reduces fast bypass in the tested packing | [Commands](adsorption/disordered_packing/README.md) |
| Adsorption | [N96 disordered packing](adsorption/disordered_refinement/REPORT.md) | The direction of the geometry effect holds and the 2-point disordered grid screen remains unmet | [Commands](adsorption/disordered_refinement/README.md) |

Filtration hydraulics is closed for now; adsorption remains research in
progress. Earlier BGK numbers and later TRT results retain separate identities.
The reports give their own experimental/scenario and operating-range limits.

Shared numerical models remain in `julia/`; shared Python runtime and interface
code remains in `src/pore2process_demo/`. Existing run folders, literature
extractions and tool installations retain their locations. See the
[layout guide](../docs/architecture/project_layout.md) and
[ordered native/process reproduction](../docs/guides/reproduction.md).

## Metadata and future examples

Every listed case has `example.toml` with its actual entry points, prerequisites,
protocol/reference files and separate validation status. Use the
[contribution guide](../docs/guides/adding_examples.md) and
[template](../templates/example/README.md) for new cases. Shared documentation is
indexed in [docs/](../docs/README.md). `scripts/check_repository.py` validates the
catalogue without executing the listed simulations.
