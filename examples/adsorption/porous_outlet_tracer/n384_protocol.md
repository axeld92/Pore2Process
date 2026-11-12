# N384 short-outlet confirmation

Use 4×4×4 tracer children of each archived N96 liquid voxel. Retain the
coarse solid mask, archived TRT flow, rectangular face reconstruction,
open axial seam and transverse periodicity. Preserve the original
4 ms pulse, 1e-9 m²/s diffusivity, 80 µs step and 10,855-step window.
The benchmark uses one porous cell.

Use the existing donor-cell and centered corrections and SSP RK3
stages. Store the reconstructed flow per original voxel and normal
subplane to reduce repeated fine-grid data access. Check that this
implementation reproduces the archived face graph and matches the
previous N288 outlet solver over short prefixes from zero tracer,
a smooth varying field and the saved late-tail field. Accumulate each
cell's face contributions in the same order.

The previous cap was 0.8 divided by maximum outgoing advective plus
diffusive depletion rate. Retain the user-specified 80 µs if its
depletion number is below 1, the convex low-order forward-Euler
limit. Check this on every fine cell and retain the original
concentration and mass screens over the full N384 history.

Compare N288 and N384 over the matched physical window. Report
time to 99% recovery, the residual outlet tail and its contributions
to moments, recovery, early arrival, mean and variance. Use the
same 99.9% recovery, last-transit 1% raw-moment, 5% mean and 10%
variance screens. Relative mesh changes use N384 in the denominator.
No additional hydraulic or full-bed execution is included.
