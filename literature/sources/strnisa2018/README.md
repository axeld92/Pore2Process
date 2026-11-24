# Strniša et al. (2018): scale-up measurements for a micro packed bed

[Chem. Eng. J. 350, 541–550](https://doi.org/10.1016/j.cej.2018.05.028)
· [metadata](metadata.json) · [evidence](evidence.toml).
The paper and its two-page supplement are both in the private inbox.

## Source separation

Three experimental datasets use the same sieved commercial catalyst
(reported mean diameter 372.10 ± 27.05 µm, p. 545), but not the same
fluid or measured quantity:

- **Hydrodynamics, tracer (pp. 543–546, Table 1, Fig. 3):** 10 µL pulse
  of 50 mM glucose in 100 mM PBS at 24 °C, flow 100 µL/min. The glucose
  detector has finite chip volume (4 µL), and Fig. 3b includes the
  response of the measurement loop with no reactor and with reactors A
  and C. All three curves are needed to compare measured transit-time
  distributions without confusing tube/detector response with packed-bed
  dispersion. The reported Table 1 moments come from the pulse experiments;
  the reported dispersion coefficients are *derived* from those moments,
  not separately measured molecular diffusivities. The full Fig. 3b
  marker curves have now been digitized in [rtd_fig3b.toml](rtd_fig3b.toml);
  Fig. 3a's peak-normalized signal cannot establish mass recovery.
- **Hydraulics (p. 544, p. 547, Fig. 6):** water at 24 °C, not the
  reacting organic mixture. The sensor specification and reported maximum
  are transcribed. The near-overlaid pressure curves merit a measurement-
  chain resistance audit before assigning all pressure loss to the beads.
  An equivalent n-heptane pressure drop was estimated by the authors,
  not measured. No blank pressure curve is reported.
- **Reaction (Table 2, Fig. 5, pp. 547–548):** vinyl butyrate and 1-butanol,
  600 mM equimolar feed in n-heptane at 24 °C. This transesterification
  produces butyl butyrate and acetaldehyde, unlike pentyl propanoate's
  water-producing reaction. No time-on-stream residual activity series is
  provided. The batch initial-rate assay (pp. 543–544) measures an
  independent activity 13.31 ± 0.27 U/mg immobilized preparation in the
  same reaction at 24 °C, 800 rpm and 5 g/L catalyst. It is one feed-rate
  anchor, not a complete rate law over substrate/product concentrations.

## Inventory and conversion

Table 2 gives empty and inferred void volumes. Recover the approximate
weighed preparation mass by inverting the authors' void-volume convention,
`W = 0.555 mg/µL × (V - Vv)`. The masses are reconstructed from rounded
values, not extra independent measurements. Check the resulting loading
against the table's U/µL **void**-volume column using the independent
specific activity. The assumed 10 wt% protein fraction is not another
mass measurement. A–C are basic channels; D–F have triangular distributors
and D/F change depth. The widths in those rows refer to the rectangular
part of the channel.

Selected markers from the three panels of Fig. 5 (native Im0; 1008 × 2367
pixels) retain the panel, symbol, y-axis calibration and nominal residence
time. Nominal time is `Vv/Q`, not catalyst age. A is the reference curve;
B and C are taken from (a), E from (b), D and F from (c). D/B recur in
other panels and are not counted as independent repeats; a few obscured
points are omitted. Near-overlapping high-conversion symbols are approximate.
1.5 percentage points is a marker-reading scenario, not experimental
uncertainty. The paper states at least two successive outlet samples at
steady state, but does not furnish numerical uncertainty for every point.

Table 2's BPN at 1.4 min is based on butyl-butyrate concentration per
mass of CALB in the void-volume convention. Reconstruct the independently
reported production rate using BPN × CALB mass / nominal time. Keep it
separate from the stoichiometric product flow inferred from the
butanol-conversion markers and feed rate; they are not interchangeable
measurements.

## Geometry for a potential next benchmark

The main article's historical numerical experiment used custom D3Q19 BGK
flow and modified D3Q7 BGK tracer, 2048 × 256 × 32 nodes, 32-node spheres
and a 6 s injection (pp. 544–545). The authors reported roughly 1% mass
imbalance and local solid-wall concentration accumulation and proposed
improving the collision scheme/walls (p. 546). Those are prior results, not
our verified solutions.

Supplement Fig. S1 shows four synthetic realizations, with porosities
0.570, 0.578, 0.573 and 0.572. Fig. S2 is a lower-resolution visualization
with 1024 × 128 × 16 nodes and 16-node spheres, not a grid-convergence
dataset or experimental specimen geometry. Neither supplies the actual
3D positions of beads in the reactor. The measured channel depth is
400 µm versus mean bead diameter 372.10 µm; the historical mesh put one
sphere diameter across the depth. Adhesion, wet size and exact clearance
are not measured, so assess a labelled clearance/size/porosity sensitivity
instead of treating those structures as exact reconstructions.

A later implementation should use the existing stock OpenLB TRT flow
and independently conservative Julia tracer transport, conserve absolute
tracer mass, and qualify the same blank-loop treatment on both channel
lengths before attempting reacting-liquid predictions. The
[microreactor tracer case](../../../examples/biocatalysis/strnisa_tracer)
implements that separate comparison. The related same-reaction
[Pohar et al. (2012)](https://doi.org/10.1016/j.cej.2012.02.035)
paper is a candidate for kinetics; its principal ionic-liquid conditions
must be checked before importing parameters into the n-heptane case.

## Complete loop-response extraction (Fig. 3b)

Native Im0 (PDF p. 6/article p. 546, 1019 × 1598 px) contains the
blank square, short-reactor triangle, and long-reactor diamond pulse
curves. The axis centers are x=107 at 0 min and x=909 at 6 min;
y=1503 at 0 mM and y=939 at 4 mM. Marker centers are approximately
11.14 px apart (5 s). The saved readings use all distinct plotted
markers from the pre-pulse region through the end of the plot
(6.67 min), and preserve the overlapping tail readings separately
from the source's uncertainty. Index 31 has no resolved markers
in the raster. A colored marker-center overlay was checked against
the source image. Interpolation between measured points is piecewise
linear; the analysis uses an added zero anchor at injection time.

The full-loop signal incorporates the 10 µL injection, tubing,
detector and packed channel. Under an unchanged linear measurement
chain, subtract blank mean/variance from the full measured moments
to infer the incremental reactor contribution. This subtraction
is a conditional measurement interpretation, not a direct
pore-dispersion measurement. The raw plot has nonzero baseline
and overlapping tail markers at the pixel scale; baseline offsets
of 2 or 5 px are separate reading sensitivities, not measured
noise statistics. Fig. 3b inset filled-marker error bars read
approximately 0.17–0.25 min around A's 0.21 min and 0.45–0.61
min around C's 0.53 min; the caption does not specify their exact
numeric values or statistical definition. Preserve the printed
Table 1 moments and these graphical ranges separately from the
Fig. 3b digitized moment differences.

### Cross-check of the dimensionless comparison (Fig. 3a)

[Selected Fig. 3a marker centers](rtd_fig3a_audit.toml) correspond to
the measured A intensity sequence in Fig. 3b. The 21 centers imply a
time scale of 118.91 s and peak normalization within 0.0033 of the
Fig. 3b values. The integrated mean of the **complete measured A loop**
is 119.04 s, while Table 1 reports 12.6 s for A. This is a reading-based
inference about the plotted curve, not access to the authors' data
processing. A normalized overlay of the experimental and numerical
curves cannot on its own establish the absolute reactor-only residence
time or dispersion. Use the dimensional blank, A and C curves for this
qualification.

For a specific subsequent apparatus check, see the [related 2018 detector
and additional tracer evidence](related_detector_evidence.md) from
Panjan's doctoral thesis. It is not a new fit to the held-out C curve.

The article's assembly details are separately retained in
[plumbing_geometry.toml](plumbing_geometry.toml). They independently
confirm Bajić's 0.50 mm tube ID; neither paper reports wet port bore
dimensions or relative bypass/bed tube lengths. The
[plumbing audit](../../../examples/biocatalysis/bajic_plumbing_audit/REPORT.md)
tests the straight-tube length hypothesis without changing the earlier
TRT or tracer predictions.
