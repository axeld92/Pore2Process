# Exact population balances

`population_audit.h` observes the installed OpenLB populations; it implements no
collision, boundary condition or transport update. It is restricted to the
single-block triply periodic cases used in this investigation. Material 2 is
the full-way BounceBack wall; all other materials are fluid.

Call `begin(block, geometry)` before stock collision, `afterCollision(block,
geometry)` after the ordinary PostCollide communication, and `afterStream(block,
geometry, step, csv_path)` after stock streaming. Root application code chooses
the recorded steps. The diagnostic does not change a population or velocity.

For each node, the expected next density is the sum of the actual post-collision
populations at its upstream neighbours. The signed flux on an outgoing link is
its post-collision population minus the opposing population at the destination.
Both definitions use OpenLB's descriptor velocities/opposites. The constant
population shift cancels in link differences; positivity diagnostics explicitly
restore OpenLB's quadrature weights.

The CSV density columns retain the scalar-oriented `concentration_*` spelling.
Their values are the native `lbm<D>::computeRho` moment: kg/m3 concentration for
the scalar, dimensionless lattice density for water. For water, analysis uses
`density_scale=1000` kg/m3. `net_to_*_per_s` is native density times m3/s; multiply
by that density scale for kg/s. Momentum and individual populations remain in
lattice units. Indexing is zero-based and coordinates are periodic.

`analysis.summarize` checks local collision conservation, exact streaming and
link balance, cancellation across fluid-wall interfaces, and total periodic
link cancellation. It distinguishes wall-adjacent fluid, bulk fluid and wall
nodes. The wall's stored populations must be included: instantaneous transfer
from physical liquid into full-way BounceBack storage need not vanish, although
the fluid-plus-wall discrete inventory remains conserved.

These balances answer whether numerical mass is being created or lost. They do
not establish bounded concentrations, correct spatial transport, or equality
between the carrier's population flux and a separately evolved scalar flux.
In particular, a centred finite-volume divergence of exported velocity is a
different diagnostic and is never substituted for these population balances.

`uniform_advection_audit.h` is a second read-only diagnostic, using stock
`equilibrium<D>::firstOrder` with a hypothetical initially uniform concentration
of 1. It sums the equilibrium populations arriving from actual upstream cells,
using their initialized `VELOCITY` field. With zero velocity in solid cells this
predicts exactly the first streamed step of a uniform scalar because both bulk
collision and wall BounceBack leave its initial equilibrium unchanged. No
velocity gradient or stencil is recreated outside OpenLB. The resulting
`scalar_equilibrium_divergence_per_s` is a scalar equilibrium-stencil property,
not the water population balance. `summarize_uniform_stencil` can compare the
prediction with the first step of an independently executed constant control.
