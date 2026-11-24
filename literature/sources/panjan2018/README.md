# Panjan (2018): dimensional microreactor tracer and detector evidence

[Peter Panjan, doctoral dissertation, Acta Universitatis Ouluensis D 1453](https://oulurepo.oulu.fi/bitstream/handle/10024/36333/isbn978-952-62-1819-9.pdf?isAllowed=y&sequence=1), ISBN 978-952-62-1819-9 (PDF). The user supplied the original PDF in the private inbox. The file SHA-256 is recorded in [metadata](metadata.json); the PDF itself is not distributed in Git.

## What is extracted

- [Conditions and detector data](evidence.toml): Sections 4.1.1, 4.3.1,
  4.4.1, Figure 23 (apparatus), Figure 25 (pulse results), and the separate
  glucose detector's Figure 4/Table 1. The PDF page number is two greater
  than the printed page number. The 4 µL detector chamber is not the
  same as a dynamically measured detector time constant. The reported
  *leading-edge* residence delays are not first moments.
- [Figure 25 trace coordinates](rtd_fig25.toml): native image Im30,
  1442 × 846 pixels, printed p. 78/PDF p. 80. Circle, square, triangle,
  diamond and line centres were read from the original image. Time
  markers are 5 s apart; the x spacing inferred in the independent
  circle/square peaks is about 15.13 pixels per 5 s. Concentration scale
  is 5 mM per 703 pixels. Overlap indices and deliberately obscured
  markers are listed; connecting-line readings are approximate. The
  extracted centres were checked over the source image. No unprovided
  replicate traces or measurement confidence intervals are inferred.
- The visible blank and R1 records finish near zero at approximately
  250 and 300 s; R2 is plotted to 400 s. The wide R3/R4 signals remain
  appreciably above zero at the 420 s plot boundary. Treat those as
  incomplete pulse integrals, not complete RTDs. R3 starts above the
  baseline before the main rise; its source and experimental sequence
  are not specified. The reported volumes for the wide configurations
  do not by themselves establish per-bead excluded volumes.

## Source independence and limitations

At clear matched raw markers, the thesis and Strniša article blank
curves agree to 0.0071 mM RMS; R1 and article A agree to 0.0133 mM RMS.
Both are within the combined native-figure pixel-reading allowance
without adjusting time or concentration. This strongly suggests shared
traces, not independent repeats. Some late-marker and window differences
are retained; do not merge their separate moment reports or average them
as experimental replicates. R2's markers differ from article C and give
an additional length comparison, but their proximity in dimensional
mean arrival is itself an important audit result.

The separate detector characterization reports a 58 s step T90 and
static calibration saturation at higher glucose concentrations. The
text on p. 72 indicates steps were taken to stay within the sensor range,
but the actual electrode calibration and complete step response for
each Fig. 25 pulse are not printed. Do not correct the pulse by assigning
the earlier calibration to it or by adding another sensor convolution
to the already measured blank.

The thesis references [Bajić et al. (2017)](https://doi.org/10.1016/j.procbio.2016.09.021)
as a reactor-design precedent; that paper uses a different immobilized
preparation. It is a potential source of separately dimensioned fittings,
not a justified substitute for the present bead geometry or volume.
See the [completed dimensional audit](../../../examples/biocatalysis/panjan_tracer_audit/REPORT.md)
for all window/reading checks, the conserved-inventory test, and the
separate reduced predictions derived from saved stock TRT/Julia ages.
