# OpenLB code boundary

The project requirement is to reuse the existing OpenLB implementation and
follow its example structure, with every project addition outside the original
distribution. This applies to future cases as well as the hydraulic benchmark.
It is also recorded in the repository's `AGENTS.md`.

As of 9 November 2026, porous-media cases must use **at least TRT** collision
dynamics. Current pressure-driven cases use stock TRT; the star channel combines
stock TRT collision and stock exact-difference forcing through OpenLB's `Tuple`.
No new collision or forcing implementation is introduced. Older BGK numerical
references remain historical and are not relabelled as TRT runs. The
[new showcase](star_channel_showcase.md) records the fresh verification.

The reference structure is OpenLB 1.9r0's
`examples/laminar/poiseuille3d/{poiseuille3d.cpp,case.h}`. Our external case uses
`createMesh`, `prepareGeometry`, `prepareLattice`, `setInitialValues`, `simulate`
and `getResults`, called in that order from a small main program. The physical
problem uses two prescribed pressure planes, rather than the example's
prescribed inlet velocity; this is an application-level setup choice.

| Responsibility | Existing OpenLB component |
| --- | --- |
| Case, mesh, geometry, lattice and vectors | `Case`, `Mesh`, `SuperGeometry`, `SuperLattice`, `Vector` |
| Velocity set and collision | Stock D3Q19, `TRTdynamics`, or a suitable more advanced stock model; stock `Tuple` composition where needed |
| Geometry | `IndicatorCylinder3D`, `IndicatorCuboid3D`, `IndicatorCircle2D` |
| Units | `UnitConverterFromResolutionAndRelaxationTime` |
| Walls and pressure planes | `setBouzidiBoundary`, `BouzidiPostProcessor`, `boundary::LocalPressure` |
| Initialization | `AnalyticalLinear3D`, `AnalyticalConst3D`, `momenta::setPressureAndVelocity` |
| Field access | `SuperLatticePhysVelocity3D`, `SuperLatticePhysPressure3D`, `SuperLatticeDensity3D` |
| Flux/pressure diagnostics | `Hyperplane3D`, `HyperplaneLattice3D`, `SuperPlaneIntegralFluxVelocity3D`, `SuperPlaneIntegralFluxPressure3D`, `SuperPlaneIntegralFluxMass3D` |
| Steady-state monitor | Two `util::ValueTracer` instances for midplane flow and tap pressure drop |

There are no project dynamics, descriptors, custom indicators, templated
numerical classes, distance solvers, quadrature implementations or boundary
post-processors. The small configuration and diagnostic records, coordination of
the stock convergence monitors, and SI CSV/JSON serialization belong to the application.
The previous `ExactCylinder` implementation was removed. Its historical results
are preserved in Git and [the earlier report](hydraulic_benchmark_legacy.md).

Half-cell transverse grid centers avoid placing an even-resolution grid point
exactly on the circular wall. OpenLB's stock cylinder distance method is called
unchanged. OpenLB's plane-integral implementation checks interpolation neighbors
even in discrete mode. We supply the surrounding material set and restrict the
integration plane with its existing circle indicator, so those checks do not
trim valid fluid-center samples. There is no analytical area correction.

All tracked application files are in `hydraulic/openlb/`. The external build
configuration is there too. Sources/headers are copied into
`.local/hydraulic/build/`; core objects/archive and case objects/executable stay
there. The five upstream core translation units are compiled directly without
modification. The case uses the upstream `onlysample` target, which avoids
in-tree core/dependency build targets. No file is written to the OpenLB tree.
The earlier 2D compatibility probe now also builds byte-identical copies of the
stock source in `.local/compatibility/openlb_example_build/`.

Before and after builds, original distribution files are compared with the
pinned official archive. The existing vendor `config.mk` from the earlier
compatibility milestone is ignored and never loaded or rewritten. Original
numerical/library/example source and build-rule files must match the archive.
New build flags, configuration and version labels are supplied externally.
Compiler/dependency paths in this configuration are specific to the verified
Apple Silicon Mac and need deliberate adaptation on another platform.

Official background: [OpenLB user guide](https://www.openlb.net/user-guide/).
The installed 1.9r0 source and examples supply the exact API reference; the
public guide page currently lists earlier releases.
