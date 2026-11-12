# Second-order-equilibrium candidate: tested and rejected as a correction

Five actual N32 OpenLB runs tested a composition of stock second-order
equilibrium, TRT, external velocity and CombinedRLB. **It does not fix the
observed porous or smooth-shear concentration errors.** Physical diffusivity,
time step, mesh, supplied velocity and observation window were unchanged.

| Actual candidate case | Minimum C | Maximum C | Final profile RMS / amplitude |
| --- | ---: | ---: | ---: |
| Porous, sine profile | .118009 | .882681 | No exact porous profile available |
| Porous, initially constant .5 | .497217 | .501831 | Constant preservation fails |
| Smooth shear, sine profile | .202785 | .797215 | 38.586% |
| Smooth shear, initially constant .5 | .500000 | .500000 | Constant preserved within 1.3e-13 kg/m3 |
| Uniform velocity, sine profile | .250541 | .749459 | .22488% |

Concentrations are kg/m3; ranges cover the whole native time history. All five
cases conserve the recorded inventory within 9.6e-13 relative. The smooth-shear
and uniform-velocity comparisons use the same independently converged continuum
reference as the unchanged parent study; uniform flow additionally agrees with
the analytical advection/diffusion expression. The uniform test passes the 1%
profile criterion, showing that the candidate is wired correctly in this
simple case. It does not qualify the porous or shear application.

The parent N32 porous range is .117696–.882507 and the parent shear profile RMS
is about 38.58% of amplitude. The candidate's effect on those failures is
negligible. Its porous constant error is below the loose engineering allowance
but fails exact constant preservation. A matched whole-window parent constant
comparison is needed before claiming a quantitative improvement in that case;
the initial two-step parent control only tests the exact first-step stencil.

Only existing OpenLB facilities are composed; original OpenLB files are
unchanged. Initial populations use stock `equilibrium<D>::secondOrder`
explicitly because CombinedRLB's ordinary initialization is hardcoded to first
order. The FirstOrder uniform-stencil header is retained as an unused copy from
the parent application; this candidate does **not** call it or produce a
`uniform_stencil.csv`, and no FirstOrder stencil prediction is used to assess
the SecondOrder runs. The legacy `initial_profile` label in porous `raw.json`
has a fixed .25 amplitude spelling; actual amplitude is preserved in
`request.json` and `peak.json` and the analysis uses that input snapshot.

The successful native cohort is `20261110T210039Z-4ac42cbc`. Its checker validates
32 source files and 363 artifacts, reanalyzes the five completed cases, verifies
all 2,070 original distribution files, and launches zero new simulations.
One preceding build attempt failed on a pointer overload in application
initialization; passing OpenLB's existing Vector fixed it. That attempt produced
zero native simulations and its source snapshot/build log is retained.

From the project root, with the parent's archived N32 request present:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.flux_audit.candidate.run
uv run python -m examples.adsorption.transport_root_cause.flux_audit.candidate.run --check
```

The runner records all actual native commands and explicit source/field paths.
It deliberately uses the preserved N32 parent request identified in `run.py`;
it does not silently substitute a velocity source on a fresh machine. No
adsorption solve, DWSIM execution, coefficient calibration or accepted transport
correction results from this candidate. See the [prospective protocol](protocol.md)
and [local population findings](../findings.md).
