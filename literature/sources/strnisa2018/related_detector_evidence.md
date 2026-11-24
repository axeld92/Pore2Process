# Related published measurement-loop evidence

**Not used to calibrate the frozen A-to-C prediction.** This source was found
while the held-out stock TRT carrier was running. It gives a concrete next
test of the apparatus and geometry, rather than a new fitted transport law.

[Peter Panjan, *Innovative microbioreactors and microfluidic integrated
biosensors for biopharmaceutical process control*, Acta Universitatis
Ouluensis D 1453 (2018), ISBN 978-952-62-1819-9 (PDF)](https://oulurepo.oulu.fi/bitstream/handle/10024/36333/isbn978-952-62-1819-9.pdf?isAllowed=y&sequence=1).
Panjan is a coauthor of Strniša et al. (2018).

## Detector characterization (pp. 42–43, Table 1)

The thesis reports a step-response $T_{90}$ of 58 s for a characterized
electrochemical glucose detector, with a 0–3 mM linear range in the table,
0.006 mM detection limit, and intra- and inter-sensor repeatability
variations of 4.0% and 5.7%. The article's stated linear range extends to
5 mM; do not assume the thesis detector's calibration/dynamic response is
identical to every article pulse. The physical detector chamber in the
thesis experiment has a reported volume of 4 µL (p. 78).

An illustrative first-order instrument with a 58 s $T_{90}$ would have a
25.2 s time constant and thus a 635 s² variance. That is comparable to the
complete blank variance of the article (769 s²). It is a scale check only;
neither a response kernel measured on the article's actual detector nor
a reason to add a second detector convolution to the measured blank.

## Additional tracer experiments (pp. 75–78, Figs. 23 and 25)

The thesis reports an aqueous 0.1 M PBS test with 100 µL/min flow and
a 10 µL, 50 mM glucose pulse. The no-reactor control connected the inlet
and outlet tubing directly. The geometries and reported volumes are:

| Thesis reactor | Length × width × depth (mm) | Reported volume (µL) | Reported onset delay (min) |
| --- | ---: | ---: | ---: |
| R1 | 25 × 3 × 0.4 | 17 | 0.21 |
| R2 | 52 × 3 × 0.4 | 35 | 0.35 |
| R3 | 52 × 16 × 0.4 | 266 | 2.61 |
| R4 | 52 × 29 × 0.4 | 583 | 4.15 |

The reported delay is based on the leading concentration rise relative
to the no-reactor control, **not a first-moment mean residence time**.
Do not mix it with Strniša's printed moments. R1 and R2 appear closely
related to the article's narrow reactors A and B. The thesis's R3/R4
volume/geometry and any overlap in trials are not established; no
additional independent replicate is assumed. The full Fig. 25 curve
would permit a useful measurement-loop transferable-length comparison,
and a separate width/distributor test if the specimens' volume definitions
and geometry are established. The online primary text was checked;
the remote page-image fetch failed and the Fig. 25 curves have not yet
been digitized. The article's Fig. 3a/3b audit remains the actual measured
comparison in the current example.
