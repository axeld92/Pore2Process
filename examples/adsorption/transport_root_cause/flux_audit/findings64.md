# N64 pore-tracer peak: local geometry and exact flux audit

The matched-time-step N64 result remains unqualified: its largest concentration
excursion is .12854336 kg/m3 outside the initial .25–.75 interval. The selected
snapshot is step 7510. Its concentration range is .12709329–.87854336 kg/m3;
the whole-history minimum, .12708898, occurs at step 7497. Location statistics
below concern the selected snapshot, not the union of all errors over time.

## Where the largest peaks sit

All 16 nodes whose excursion exceeds .1 kg/m3 are next to solid and near the
downstream poles of the declared smooth spheres. Their smooth-surface outward
normal has x component +.99927034, and the extreme nodes have nine neighbouring
solid links. The local speed is .000171546 m/s, much smaller than the mean
interstitial speed. The peak-node signed distance from the smooth sphere is
positive: **1.85091 micrometres, or .08527235 grid spacings**.

| Mesh | Selected maximum C, kg/m3 | Peak-node gap, micrometres | Gap / dx | Local speed, mm/s |
| --- | ---: | ---: | ---: | ---: |
| N32 | .88250691 | 13.5514 | .312161 | .755190 |
| N48 | .86812205 | 20.1446 | .696054 | .498936 |
| N64 | .87854336 | 1.8509 | .085272 | .171546 |

The placement of the last fluid nodes relative to the smooth surface changes
substantially with mesh. This accompanies the nonmonotone extrema and supports
a separate grid-phase sensitivity test. It is **correlation, not proof of a
cause or a monotone gap/error law**: N64's gap is smallest, but N32's excursion
is slightly larger. The full-way numerical wall is represented by voxel
midlinks, not by the analytic spherical cut; these gap values must not be
described as the exact distance to the discrete bounce-back wall.

Smaller violations extend into the liquid. N64 has 2,424 fluid nodes outside
the engineering allowance in this snapshot: 80 wall-adjacent and 2,344 with no
solid in their immediate 26-neighbour geometric neighbourhood. The physical
solid surfaces are convex spheres. The snapshot does not establish where an
error first developed or how it propagated.

## Mass conservation and population positivity are different checks

The deterministic peak replay matches every saved concentration **exactly**.
Local collision density changes are below 5.6e-16; actual population streaming
and signed-link balance residuals are below 8.9e-16 native concentration units.
Fluid-wall exchanges cancel when numerical wall storage is included. Thus
the saved concentration distribution is not explained by a mass leak or an
incorrect diagnostic reconstruction.

Unlike the N32 peak, N64 has **16 negative unshifted post-collision populations
at 16 fluid nodes**, with minimum -1.15845e-5. These nodes do **not** coincide
with the 16 severe concentration peaks. They have C≈.227220 kg/m3, 13 solid
neighbours and nearest-surface normal x≈-.94526 on the upstream side. Their
smooth-surface gap is also .08527235 dx. None of the severe downstream peaks
has a negative post-collision population in this audited step. Negative
populations are an additional warning in this N64 case, not a sufficient
explanation of the worst extrema; N32 already violated concentration bounds
without any negative populations at its selected peak.

The [machine-readable record](findings64.json) retains the native summary/
artifact hashes, precise locations, negative-population node list and exact
balance diagnostics. It does not add simulations, fit a coefficient or qualify
an adsorption/process result. Replay with:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.flux_audit.findings64 --check
```

The separately verified [flat-wall control](flat_channel/README.md) resolves its
profile below 1% RMS at N96 with unchanged stock dynamics and D. That establishes
that resolution can fix a simpler wall-bounded problem; it does not establish
how much refinement or what geometry treatment the FCC case needs.
