# Stock TRT stationary-wall comparison, 11 November 2026

The [scientific report](../../../examples/adsorption/stock_transport_candidate/REPORT.md)
records six matched native OpenLB executions against the converged Julia
radial reference from the stationary-wall benchmark. Plain TRT, half-way
reflection and their combination all failed uniform-tracer and profile
screens. The closest alternative was regularized TRT with half-way
reflection, with 11.82% final profile RMS error and a 0.06714 kg/m³ uniform
departure. Relative closed inventory drift was below 6.33e-13 in all cases.

The run `20261111T133159Z-7f69a322` passed read-only replay of 32 sources
and 168 artifacts. The [run guide](../../../examples/adsorption/stock_transport_candidate/README.md)
records commands and archived inputs. The study generated zero bed or
process simulations.
