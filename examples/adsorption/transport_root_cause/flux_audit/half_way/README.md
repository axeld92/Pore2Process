# Stock half-way scalar reflection candidate

This separate, finite diagnostic changes the scalar wall treatment while
retaining the physical velocity, diffusion, TRT collision and comparison
window. It uses OpenLB 1.9r0's existing `setBouzidiBoundary` registration and
`BouzidiPostProcessor` with active distances set to exactly 0.5. No boundary,
collision operator, descriptor or original OpenLB source is reimplemented.

The candidate keeps the voxel midlink wall convention. It does not use the
analytical sphere intersection distances for curved-wall interpolation. The
stock analytical indicators are used to register links; only already active
distance entries are configured to 0.5. Inactive entries remain -1.
The installed `src/boundary/setBouzidiBoundary.h` supplies the operator and
`examples/laminar/poiseuille3d/case.h` shows its stock registration pattern.
The OpenLB maintainers describe the full-way/half-way distinction in their
[boundary discussion](https://www.openlb.net/forum/topic/bounceback-boundary/).

The stock PostStream operator immediately returns outgoing populations to the
fluid. Solid nodes use NoDynamics; their streamed population copies are not
additional physical inventory. The first two and final steps have read-only
reflection and local fluid balance diagnostics. Link tables include every
active direction and periodic seam. Initial equilibrium is restored on fluid
after stock initialization, then checked against the declared concentration.

The [prospective protocol](protocol.md) requires four flat controls before two
porous tests. The flat comparison uses the independently qualified continuum
oracle retained in `../flat_channel/reference_qualified.py`: one Fourier mode
in the streamwise direction and Neumann diffusion across the channel. Physical
wall locations are 0 and H, fluid centers are half a cell from each wall, and
there are four periodic z layers. Earlier full-way controls explicitly tested
the equivalence of four and 32 layers for this z-invariant problem.

Run from the repository root:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.flux_audit.half_way.run
uv run python -m examples.adsorption.transport_root_cause.flux_audit.half_way.run --check
uv run python -m examples.adsorption.transport_root_cause.flux_audit.half_way.independent_review \
  --summary results/adsorption/transport_root_cause/half_way/runs/20261110T220640Z-532998c7/summary.json
uv run pytest tests/adsorption/test_transport_half_way_review.py -q
```

The runner preserves commands, process IDs, exit codes, inputs, source snapshots,
binary hashes, native fields, link tables and results in a new timestamped
directory under `results/adsorption/transport_root_cause/half_way/`. Its checker
reassesses those artifacts without running another simulation. The independent
review is separately recorded and does not mutate the native simulation inputs.

## Result: rejected as a particle-bed correction

All six native cases completed. The flat stage passed its declared gates. A
constant remained constant to 6.9e-14 kg/m³, the diffusion-only profile had
0.228% RMS error relative to its initial amplitude, and the N96 Poiseuille
profile was bounded with 0.697% RMS error. The fine global RMS result does
not imply better wall accuracy: the first fluid row had 3.59% RMS error and
the maximum pointwise error was 5.07% of amplitude. The corresponding full-way
values were 0.896% and 2.26%.

| Particle mesh | Original excursion, kg/m³ | Half-way excursion, kg/m³ | Half-way all-time concentration range, kg/m³ |
| --- | ---: | ---: | ---: |
| N32 | 0.132507 | 0.169700 | 0.080918–0.919700 |
| N64 | 0.128543 | 0.176602 | 0.080269–0.926602 |

The largest excursions grew by 28.1% and 37.4%, respectively. The native
diagnostics close the actual reflection and local fluid balance at the
recorded steps. Physical fluid inventory drift is below 7.5e-13 in all cases.
The separate review matches all reflecting links against the exported material
and descriptor, including periodic seams, and checks the actual initial fields,
fluid inventory, diffusivity, and unchanged porous velocity and run arguments.
The candidate checker replayed 34 sources and 441 artifacts. Four independent
topology tests passed. No vendor source changed.

The [peak localization](localization.json) puts the largest N64 maximum at
the same node as the full-way run. The count of fluid nodes exceeding the
physical bounds by more than 0.1 kg/m³ at the selected peak increases from
16 to 48; they are all beside downstream particle poles. This association
does not identify the first source of the error.

This rejects immediate return at the same voxel wall as a remedy for these
particle-bed peaks. It neither tests a smoothly interpolated curved wall nor
qualifies adsorption, residence-time distribution, dispersion or service time.
