# Filtration evidence: roles of the three supplied papers

Reviewed 8 November 2026. Main-article Markdown and conversion metadata for
Mo (2022), Biesheuvel and Verweij (1999), and Jan et al. (2025) have been read.
The relevant figures were inspected. Full texts and images stay in the ignored
`inbox/` folders. The [source ledger](../../sources/index.json),
[Mo extraction](../../sources/mo2022/extraction.json), and
[support/membrane extraction](support_membrane_extraction.json) retain provenance.

| Paper | What we can use | What it does not establish |
| --- | --- | --- |
| Mo et al. (2022) | Fixed-top-layer support plateau; separate intermediate resistance; model settings; supplied supplement now resolves S2/S3 citations | A fully specified matched layer/geometry dataset; original experimental sources still need audit |
| Biesheuvel and Verweij (1999) | Nine measured support-property cases, processing details, and a permeability/strength design framework | A matched measurement of support, intermediate and selective resistances in the same complete membrane |
| Jan et al. (2025) | Experimental whole-membrane permeance, coating morphology and ion-rejection context | Support-only improvement with the selective layer unchanged |

## Biesheuvel: an experimental support reference

Table 2 reports the following support properties. The original transcription
from the supplied Markdown has now been visually corroborated against the
same study in the public author thesis, Chapter 7, printed p. 105; see
[the audit](../../sources/biesheuvel1999/source_audit.json). These measurements are accepted for
a limited reduced support closure. [The completed calibration](../../../docs/verification/filtration/calibration.md)
improves error but fails its primary 25% screen. The original OpenLB geometry
and complete membrane have not been experimentally validated.

| Row | Powder | Sinter temperature (°C) | Hold (h) | k (10^-16 m²) | Porosity | Reported mean strength (MPa) | Weibull m |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | CR-6 | 1260 | 3 | 6.99 | 0.400 | 97 | 10.7 |
| 2 | AKP-30 | 1050 | 1 | 1.90 | 0.349 | 109 | 9.8 |
| 3 | AKP-30 | 1100 | 1 | 1.35 | 0.286 | 199 | 12.5 |
| 4 | AKP-30 | 1150 | 1 | 0.766 | 0.237 | 146 | 5.8 |
| 5 | AKP-30 | 1100 | 0 | 1.35 | 0.315 | 134 | 5.0 |
| 6 | AKP-30 | 1100 | 5 | 1.37 | 0.284 | 201 | 6.3 |
| 7 | AKP-30 | 1100 | 10 | 1.26 | 0.273 | 140 | 5.3 |
| 8 | AKP-15 | 1100 | 1 | 3.89 | 0.297 | 144 | 12.1 |
| 9 | AKP-50 | 1100 | 1 | 1.14 | 0.299 | 92 | 9.7 |

Section 3.2 identifies **distilled water as the permeating fluid**. Nitrogen
pressurises the water cell; the table is not a nitrogen-permeability dataset.
The cell was pressurised up to 2 MPa, with measurements at at least three
pressure differences. Actual per-specimen pressure/flow series and water-test
temperature are not supplied in this main-article extraction. The model
viscosity in Table 3 must not be silently assigned to the measurement protocol.

The reported strength comes from four-point bending and Weibull analysis,
using 15 specimens per material type (Section 3.3). The Weibull modulus is not
a standard deviation of permeability or a ready-made device pressure rating.
Specimen size, stress field and processing matter when interpreting strength.
Zero hold time in row 5 still includes the heating/cooling cycle.

The measured k range, 7.66e-17–6.99e-16 m², provides a material-specific
comparison scale. Matching one k value and porosity would not identify a
unique connected three-dimensional pore geometry. Powder diameters are not
capillary radii; agglomeration is explicitly relevant, particularly for AKP-50.

The calculations in Table 3/Figs. 3–9 are a **separate model**, maximising
allowable flow with a pressure/strength constraint. They are not the measured
Table 2 rows, nor the same objective as our fixed-flux energy question. Their
engineering lesson is useful: a more permeable or thinner support may lose
mechanical capacity. A graded support should retain a fine, coatable surface
while opening the bulk, and its feasibility needs an explicit strength check.
The paper discusses this manufacturing direction in Section 4.3.

No radial-flow or mechanical equation was added to the software. Before
reproducing that model, check the original PDF: the logarithm in converted
Eq. 12 appears inconsistent with the series expression in Eq. 13 and an
independent cylindrical Darcy derivation. Other converted units, especially
bending-test motion units, also need verification. These are source-extraction
checks, not a finding of an error in the published paper.

## Jan: a complete-membrane reference, with a changing selective layer

Section 3.2.1 reports approximately **277 LMH/bar** for the uncoated commercial
alumina UF membrane and **26 LMH/bar** after 40 minutes of SiC deposition.
The corresponding SI permeances are 7.69444e-10 and 7.22222e-11 m/(Pa s).
The article calls this quantity pure-water permeability; its reported units
identify a permeance, distinct from Darcy permeability in m².

Figure 4a confirms the comparison visually and shows error bars. Their exact
values, definition and replicate count have not been extracted. The coating
changes selective-layer pore structure and chemistry, so this pair cannot
calibrate our unchanged-selective-layer support comparison. The uncoated
commercial membrane already contains a selective layer: it is not a bare support.

Section 2.1 distinguishes tubular hydraulic specimens from flat-sheet
characterisation specimens. Tubes have 7 mm inner diameter, 10 mm outer
diameter and 10 cm length; the flat specimens have different nominal pores.
Do not combine the flat-sheet TEM coating dimensions with tubular pore sizes
to construct a supposedly measured three-dimensional geometry. SEM and mercury
intrusion likewise do not supply full three-dimensional connectivity.

The setup is constant-flux crossflow, with details referred to Jan et al.
(2024), DOI 10.3390/membranes14010022. The pure-water test temperature and
individual flux/TMP values are not confirmed here. The 20 °C mentioned for
water's dielectric constant does not establish the hydraulic test temperature.
Figures 1–5 were inspected; a separate Fig. 6 image was not present in this
export. The paper states that the authors do not have permission to share data,
so raw series should not be assumed available.

Section 3.1.1 also identifies support/selective-layer interface effects and
coating penetration as relevant to resistance. That reinforces the need to
hold the interface region fixed in our geometry comparison. This does not
quantify a standalone interface resistance for our model.

## What changes in the research plan

The software's synthetic reference remains valid. The evidence now gives it
three distinct reference roles:

1. Use **Biesheuvel Table 2** as the first support-property dataset for a
   bounded reduced-model study; the table and units are now audited. Treat
   permeability and mechanical properties as tied to its processing route.
2. Use **Mo Figure 4** for the fixed-top-layer saturation mechanism and its
   three-layer discussion to keep intermediate resistance explicit. Its model
   values are not experimental calibration. See [the detailed review](../../sources/mo2022/review.md).
3. Use **Jan** as whole-membrane/coating context and a possible later charged
   separation case, with its own specimen and test conditions.

Do not subtract a Biesheuvel support resistance from a Jan membrane resistance.
They describe different materials, specimens and preparation routes. Such a
subtraction would create an unsupported selective-layer property.

Mo's Word supplement is now fully received and reviewed; see
[the supplement evidence report](../../sources/mo2022/supplement_review.md). The priority
original experimental source was **Zhu et al. (2021), DOI 10.1016/j.ceramint.2021.01.088**.
Its supplied main article/figures now support [a bounded matched-series fit
and prediction check](../../../docs/verification/filtration/zhu_jpeg.md), with the conservative
digitization screen still failing. Mo's original PDF remains
useful for equation checks. The Biesheuvel table audit is
complete using the author thesis; its raw measurements/uncertainty and a
specimen-linked morphology would be needed for stronger validation.
Jan's supplement would clarify the interface/coating evidence and may help
with test details; its cited setup paper is a follow-up if conditions remain
unresolved. No further broad article search is needed before checking these.

With those checks complete, select one support baseline and an explicitly
identified or assumed fixed-layer resistance. Build the controlled graded
geometry using stock OpenLB components, preserving the selective-side opening,
support thickness and representative area. Report separately what is fitted,
independently validated, and only an illustrative assumption. No production
solver settings or completed simulation outputs were changed by this review.
