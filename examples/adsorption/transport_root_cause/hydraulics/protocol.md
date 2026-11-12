# Fixed-flow water refinement protocol

This study adds independently converged N64 (and optionally N96) native
OpenLB fields to the existing N32/N48 d800 synthetic FCC case. It does not
change any historical result, adsorption parameter, physical flow, viscosity,
or geometry. The parent hydraulic inputs remain the operating-point boundary.

- Particle diameter: 0.0008 m. Nominal FCC liquid fraction: 0.6.
- Physical flow: 7.3e-7 m³/s in the existing 11 mm diameter column.
- Native viscosity/density: 0.001 Pa s and 1000 kg/m³. Native velocity is
  matched to the parent carrier Reynolds number; the scalar reader applies
  the same recorded carrier/native kinematic-viscosity ratio exactly once.
- Stock D3Q27 TRT, Kupershtokh force, stock sphere indicators and BounceBack;
  geometry and lifecycle are copied from the earlier external flow case.
- Relaxation time is 0.5 + 0.0003 N, as in the earlier spatial comparison.
- First run at a new grid starts at rest with the accepted next-coarser
  pressure gradient. Subsequent pressure-root iterations use the immediately
  previous converged same-grid native field as an initial estimate only.
- Every iteration must converge natively; source velocity is never obtained
  by spatial interpolation or an assumed analytical rescaling.
- A bounded positive secant pressure iteration targets native superficial
  velocity within 1e-4 relative. Other unchanged screens are: plane mass flux
  spread ≤1e-6; final six-sample relative standard deviations ≤1.0001e-7;
  maximum Mach <0.05; density departure <1e-4; voxel fraction within 0.015 of
  0.6. The hydraulic mesh pressure difference is reported, not hidden.
- New final-step population diagnostics may observe stock collide/AndStream
  stages. These only read populations and retain their native shifted-population
  convention; they do not alter collision, force, boundaries or fields.

Outputs retain exact requests, source snapshots, binary hashes, convergence
histories, SI velocity CSV, particle/domain dimensions, realized void volume,
and parent artifact hashes. Scalar concentration convergence remains a separate
acceptance question. No DWSIM or adsorption prediction is introduced here.
