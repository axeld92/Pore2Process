# Matched constant-background diagnostic

Declared 10 November 2026 before execution. Reuse the exact archived N64
FirstOrder/CombinedRLB/TRT executable and qualified water field. Preserve the
sinusoidal parent's physical time step, duration, diffusivity and mean C=.5,
but set initial sinusoidal amplitude to zero. Capture the actual population
field at the sine parent's peak step 7510, as well as the constant run's own
peak and final fields.

The passive numerical update is linear. Compare the constant-background
departure from .5 with the original sine result at the same nodes and time.
The whole-history maximum constant departure bounds its possible contribution
to an individual sine excursion. Report the difference field only as a
diagnostic decomposition. It is not a transport correction, renormalization,
new production model or qualification of the sine case. In particular, do not
subtract it from a process calculation or change mass/flux data.

This adds one scalar execution while the independently converged N96 hydraulic
study runs. No additional grids, molecular diffusion, collision parameters,
boundary implementation or experimental fitting are introduced.
