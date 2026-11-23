# Catalyst conditioning protocol

Retain the 2025 batch-derived rate reconstruction, its apparent
time and concentration convention, the previously tested central
Wilson–Geankoplis film model, and all rate and transport parameters.
Keep all 6.0 g of weighed catalyst present. None of the 20 measured
reactor conversions may be used to set the conditioning scenarios.

## 1. Aqueous exchange timescale

Assume that a normalized blue-dextran tracer follows a mobile
liquid clock with mean `mu` and variance `v`. Sodium azide shares
this clock but also undergoes Poisson entry per unit mobile time
into an immobile domain with exponential waiting times of mean
`tau_s`. The extra mean immobile residence per unit mobile time
is `beta`. The conditional expectation and variance of the total
time are then `(1+beta)*T` and `2*beta*T*tau_s` at fixed mobile
time T. Therefore, after subtracting the same inlet-pulse mean
and variance from both tracer curves:

```text
delta_mu = mu_azide - mu_blue
beta = delta_mu / mu_blue
var_azide = (1+beta)^2 * var_blue + 2*delta_mu*tau_s
```

Apply this formula both to the published moments and to the
integrals of the digitized curves. Assume a shared rectangular
7 min inlet pulse as the nominal calculation. If the shared
pulse has support in [0, 7] min, use zero inlet mean and the
maximum 7²/4 min² variance together for a conservative upper
envelope on the inferred return time (not a physically realizable
pulse). A negative return time for any proposed pulse would
invalidate that version of this simple exchange model.
This is a conditional normalized-shape interpretation only:
absolute tracer dose closure fails and the reactive-liquid wetting
time remains unknown.

## 2. Finite-batch moisture inventory

For the independent isotherms at 40 °C:

```text
q(c) = q_max*c/(K+c)       [kg water/kg catalyst]
c                          [kg/m³ liquid]
```

At matched conversion X, calculate produced water from the
nominal initial acid moles times X. For an initially dry batch,
use the smaller of the produced inventory and
`m_catalyst*q_max` as a maximum retention *scenario*.
Evaluate a separate reduced capacity when fresh catalyst initially
contains 0.9% water (the prior independent moisture measurement's
high sensitivity). At small X, solve the closed finite-volume
equilibrium `C_total = c + m_catalyst*q(c)/V_liquid`, reporting
the fitted equilibrium only where c <= 1 kg/m³. For water
concentrations beyond the common measured range, do not report
an extrapolated sorption equilibrium.

At steady packed-bed operation, stationary sorbed inventory
does not continuously remove water from each liquid element.
For the deliberately favorable batch-retention scenario, compare
the water activities of the bed and the batch having retained water
up to the analogue saturation capacity. Use the earlier measured
liquid-activity projections and compare only inside their supported
composition range. Compare the relative differential rate of the
independent MTBE water-activity experiment only where both
activities fall within its measured activity range. Holding these
40 °C capacity analogues across the setpoints is a sensitivity,
not a measured temperature dependence.

## 3. Transient water storage within the bead

In a spherical continuum with local equilibrium, immobile sorbed
water increases the storage responding to gradients in free pore
water concentration c. With radial coordinate r and porosity eps:

```text
s(c) = eps*c + rho_cat*q(c)
ds/dt = (1/r^2) * d/dr[r^2*D_eff*dc/dr]
D_eff = D_bulk*eps/tortuosity
s      [kg water/m³ bead]
rho_cat [kg catalyst/m³ bead envelope]
```

At the bead center use zero flux. At the surface use a liquid film
flux into the bead driven by the difference between a prescribed
external liquid concentration and the bead surface concentration.
Derive the catalyst mass per envelope volume from the published
6 g, bed dimensions and aqueous interbead porosity; keep this
density fixed for both particle-porosity scenarios. The geometry
and tortuosity come from the 2025 article; the bulk water diffusivity
comes from the independently measured neat-pentanol value scaled
by T/viscosity. Use 40 °C (the isotherm temperature) and the film
at the slowest flow (0.7 mL/min), then use it for all illustrative
displacement-age averages. Treat the aqueous (0.40) and nitrogen
(0.22) particle porosities as two distinct scenarios. Count the
entire measured sorption as additional storage; overlapping
water in the pores could reduce the additional capacity.

Test steps from c = 0 to 0.42 and 1 kg/m³ and the reverse steps,
within the common concentration range of the source figures.
These model limiting low-moisture changes, not the complete
reaction concentration history. Also test a no-sorption reference.
Use a conservative radial finite-volume balance with backward
Euler time integration, Newton updates of the storage term and
zero center flux. Verify the zero-sorption, fixed-surface case
against the exact independently derived normalized mass remaining:

```text
6/pi^2 * sum_{n>=1} exp(-n^2*pi^2*D_bulk*t/(tortuosity*R^2))/n^2
```

For both wetting and drying, report the remaining fraction of
initial-to-equilibrium water inventory difference:
`(M(t)-M_eq)/(M(0)-M_eq)`. Also compute the same measure in the
outer 80 µm shell (the geometric motivation is the external
protein-shell observation in the previous [particle evidence](../pentyl_propanoate_transport/particle/evidence.json)).
These quantities describe water differences, not the fraction
of inactive catalyst or the fraction of unfilled pores. For a
bed-scale contact-age illustration, expose successive sections
of the bed as the liquid displaces one interbead void volume
and integrate the moisture curves at the corresponding ages.
This assumes a uniform final exterior concentration. It does
not model reaction-generated water or moving wetting fronts.
Require monotonicity, bounds, and closure of integrated water
exchange to 1e-7 of the step's total storage change. Require
normalized inventory differences below 0.005 at 1-minute samples
upon grid/time refinement. After the initial 32/64-cell attempt
failed this criterion, refine to 64/128 radial cells with 2/1 s
time steps, respectively.

## 4. Constant rate sensitivities and startup stability

Use the minimum of two same-progress fresh/reuse time ratios in
the 2023 acetate study for a mild apparent relative-rate case;
digitization sensitivity displaces each marker coordinate by
its recorded reading half-range. Other cases apply the residual
activity ratios calculated from the independent Meneses assays:
benzyl alcohol contact relative to fresh, the 1:5 two-cycle
case relative to fresh, and neat propanoic acid contact relative
to benzyl alcohol contact. These are controlled comparisons of
magnitude. They do not identify the true pentyl propanoate bed
activity. Apply the factor once to apparent rate (including
inside the film balance) at constant catalyst mass. Do not
combine independently hypothesized reductions to tune to the
reactor. Seal the prediction table before reading the measured
20-run conversion file. Confirm unscaled replay against the
previous film model, axial grid refinement (512/1024 cells),
acid conservation and film balance.

Separately illustrate two simplified catalyst histories at a
magnitude taken from the previous geometric rate-basis audit:
`f(30 min)=0.64`, not a fitted activity. Under a first-order loss
from fresh catalyst, the relative activity at 180 min would be
`0.64^6`. Under first-order recovery from an initially completely
unavailable state, it would be `1-(1-0.64)^6`. Compare computed
steady conversion snapshots at these factors to the 30–180 min
stability observed in that paper's Figure S2. This is a
quasi-steady scale check on possible drift, not a transient
reaction solve or a fit to the time curve. The experiment's
unknown flow and temperature prevent a condition-by-condition
comparison. Rapid irreversible loss before steady sampling, or
a changed steady intrabead environment, is not excluded.
