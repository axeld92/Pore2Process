# Can sample handling or bead diffusion explain the conversion gap?

**Bottom line.** Published liquid data give small changes in the reported
conversion. The particle scenarios anchored to the available porosity,
bead size and reference effectiveness narrow the gap at the full
6.0 g loading. Even combining the strongest anchored bead change with
complete water omission leaves a mean overprediction of 5.60 percentage points for the
2025 reactor reconstruction. The actual reacting-liquid diffusivity and
the cause of reduced bed activity remain unresolved.

## Liquid readout

The packed-bed study measures remaining acid by titration and reports
conversion from the outlet concentration relative to the inlet. Using
published pure-liquid densities and binary excess volumes at 25 °C
changes the predicted whole-liquid readout by at most 0.06 percentage
points. Omitting all produced water lowers it by at most 1.54 points.
Measured acid–pentanol–water tie lines, with the ester assigned to
the organic phase, give a pentanol-rich readout shift of up to 0.27
points over available cases. The actual phase handling, sample
temperature and four-component outlet state remain unreported.
[Liquid data and comparisons](phase/README.md)

At 60 °C and 3 mL/min, matching 49.7% reported conversion from the
60.23% molar prediction by volume change alone would require about
21% less sample volume than the inlet basis, with all residual acid
retained. This gives the scale of the discrepancy. The applied volume
and partitioning scenarios concern how conversion is read out.
Water-dependent catalytic activity is a separate physical effect.

## Internal transport

The dry and aqueous porosity estimates provide a diffusivity-ratio
scenario of 0.55 at fixed tortuosity and molecular diffusivity. The
largest sieved diameter contrast measured on the same support in
another reacting liquid adds a diameter ratio of 0.828. Together
they reduce the apparent rate to 87–89% of the reference with the
published first-order effectiveness anchor. Reconstructing the dry void
fraction from the pore volume with the reported polymer density on a
skeletal basis gives 0.182 and a relative rate of 83–86%. Using an exterior protein
shell reduces the internal penalty. All calculations retain the
weighed catalyst inventory and the prior rate calibration.
[Particle evidence and model](particle/README.md)

| Case, using the 2025 reconstruction | Conversion at 60 °C and 3 mL/min | Error across 20 runs (RMSE, pp) |
| --- | ---: | ---: |
| Full 6 g, reference | 60.23% | 10.46 |
| Printed porosity and diameter contrasts | 56.26% | 7.83 |
| Skeletal-density porosity and diameter contrasts | 54.98% | 6.99 |
| Skeletal-density porosity, diameter contrast and water omission | 53.44% | 6.08 |
| Uniform-sphere strong-diffusion limit at printed porosity contrast | 48.28% | 3.81 |
| Measured | 49.70% | |

The strong-diffusion limit shows that sufficiently severe diffusion
could matter. It assumes a reference much more diffusion-limited
than the Table B1 estimate. With that estimate and fixed diameter,
a 0.64 relative rate from bead diffusion would require effective
diffusivity to fall to 13–16% of its reference value. That target
comes from the earlier catalyst-volume comparison. The actual magnitude
of the batch-to-bed change remains unknown.
The 2026 batch-law branch also retains a gap, with 8.16 points RMSE
in the combined scenario. All 20 reactor observations were reserved
for scoring. [Fixed-inventory results](particle/reference/bed_comparison.csv)

## Next steps

Recheck liquid-film transport around the beads with a dimensionless
Schmidt number and independently sourced liquid transport properties.
The earlier paper's printed Schmidt expression is dimensionally
inconsistent. This is the next independently testable resistance
between the flowing liquid and the catalyst. If it is small, use
independent batch activity data to check how local liquid composition
and water activity change the rate. Keep the reactor measurements
held out. [Published transport equation audit](../../../literature/sources/lorenzo2025/author_model.md#mass-transfer-treatment)
