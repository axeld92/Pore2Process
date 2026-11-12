# Prospective periodic smooth-profile isolation

10 November 2026, following the unphysical open-boundary front trials. Test a
smooth initial field C=.5+.25*sin(2*pi*x/L) kg/m3 in a fully periodic domain.
No inlet, outlet, adsorption or source is present. The application samples this
analytical initial condition; all collision, geometry and wall algorithms remain
stock OpenLB. Use the same combined regularized TRT D3Q27, N32, dt, physical
velocity and Dm=1e-9 as the earlier engineering trials. Execute one free-flow
and one FCC-bed control for EACH particle-size scale for two physical unit-cell
transits. These are short isolation controls, not column or dispersion validation.

In uniform free flow the independent exact solution is
C=.5+.25*exp(-Dm*(2*pi/L)^2*t)*sin(2*pi*(x-v*t)/L).
Require final RMS error / .25 <=1% and kinetic mass residual <=1e-8. For both
domains record min/max over all time steps and require C>=.23 and C<=.77: a
chosen .02 kg/m3 allowance around the continuum .25..75 bounds. Report liquid
and numerical wall inventory separately; wall populations are not physical
liquid. A passed short periodic control cannot waive the failed fronts or supply
a dispersion coefficient. A failed porous control with passing free flow would
localize the problem to some part of the porous velocity/wall configuration,
without uniquely identifying its cause. A passing porous control would show
that smooth positive initial conditions are better behaved than tested fronts.
