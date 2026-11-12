# N288 tracer refinement at the matched N192 time step

The [protocol](n288_protocol.md) uses the same saved stock TRT field and
N96 solid voxels as the earlier [fixed-flow refinement](REPORT.md).
Each original liquid voxel contributes 27 equal N288 tracer cells,
for 14,338,080 liquid cells and 42,753,312 shared faces. The tracer cell
width is 4.824 µm. Reconstruction holds the carrier normal velocity
constant over each coarse face and interpolates along the normal inside
each voxel. The nine subfaces of every original interface recovered
its saved flow to roundoff. The maximum parent-face sum error was
2.56e-16 of the largest parent flow. The maximum plane throughput error
was 1.80e-13 relative to the archived axial flow, and liquid volume
stayed at 1.6091e-9 m³.

Both grids use 1,314 steps in each of the two archived observation
intervals, with steps of 81.7143 and 81.7287 µs. The N288 stability cap
is 107.516 µs. Molecular diffusivity remains 1e-9 m²/s and the initial
sine is sampled at each mesh's tracer cell centres.
[N192 half-step comparison](N192_TIME.md)

The threaded Julia kernel evaluates the same bounded face increments and
accumulates each cell's contributions in the original serial order.
Four matched steps from the archived N192 initial, midpoint and final
fields reproduced every concentration bit for bit. The largest relative
difference in the dissipation sums was 8.76e-12.

## Matched-time results

Changes are relative to the N192 half-step run, with the phase reported
as the shortest signed angle between fitted longitudinal waves.

| Measure | Midpoint change | Final change | Limit |
| --- | ---: | ---: | ---: |
| Liquid concentration variance | +1.61% | +2.56% | 5% |
| Fundamental amplitude | +0.082% | −1.67% | 5% |
| Fundamental phase | +0.02177 rad | +0.05584 rad | 0.05 rad |

Final variance rose from 0.01300247 to 0.01333562 kg²/m⁶. The fitted phase
rose from 2.24644 to 2.30228 rad. Variance and amplitude pass at both
observations. Final phase still misses the criterion by 0.00584 rad,
about 0.33°. Its shift is smaller than the earlier N96–N192 shift of
0.09763 rad. The prior N96–N192 comparison used a different time grid,
so the two comparisons are reported as sensitivity checks.

The estimated numerical share of final half-variance loss fell from
5.42% to 2.00%, with absolute numerical loss decreasing from 5.929e-4 to
2.158e-4 kg²/m⁶. The N288 solution stayed in 0.2500149–0.7499851 kg/m³,
with relative liquid inventory drift of 3.34e-16. The reconstructed
hydraulic field and voxel wall still come from N96. The study measures
spatial sensitivity of tracer transport on that fixed flow.

## Verification record

The N288 run completed all 2,628 steps in cohort
`20261111T164622Z-a501da76`. All stored times match the archived N192
half-step history exactly. Read-only replay checked 26 source hashes
and 40 artifacts, including rebuilding and checking the refined
face arrays. The five existing bounded-solver analytical checks passed.
This study generated one new N288 Julia tracer run, with zero new
OpenLB, full-bed, or DWSIM runs.
[Run and replay commands](README.md#n288-matched-step-spatial-check)
