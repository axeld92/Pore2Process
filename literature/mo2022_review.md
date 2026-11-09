# Mo et al. (2022): evidence review for the filtration case

Reviewed 8 November 2026. DOI: **10.1016/j.watres.2022.118658**.
The supplied main-article Markdown, conversion metadata and all seven scientific
figures were inspected. Source files remain unchanged in
`inbox/Mo2022Dissecting/`; full text and figures remain outside Git. The
[extraction record](mo2022_extraction.json) preserves file hashes and provenance.

The paper supports the resistance-partition argument behind our research
question. It does **not yet supply an accepted experimental calibration**.
Its support/top-layer matching curves are model predictions. The supplement
has now arrived as a Word file and is fully reviewed in
[the supplement evidence report](mo2022_supplement_review.md). It resolves the
experimental reference identities; it does not supply a complete matched-layer
or geometry-calibration dataset by itself.

## What the paper establishes for our next step

Section 2.1 treats bottom, intermediate and top layers as separate hydraulic
resistances in series. This agrees with our existing Julia boundary
`delta_p = mu*J*(Rs + Rselective + Rfixed)`. In a three-layer case,
`Rfixed` must explicitly include the unchanged intermediate-layer contribution.
It must not be relabelled as selective-layer resistance.

Figure 4 is the closest comparison to our question. It holds the top layer
fixed while changing the support. Permeance approaches a plateau as support
resistance becomes small. At the same fixed filtration flux, this translates
to diminishing savings in pressure and pump shaft energy. This interpretation
is our scale-translation inference, not a DWSIM result reported by Mo et al.

The plotted hydraulic result is a plateau, rather than a falling branch that
identifies a unique unconstrained hydraulic optimum. The useful engineering
choice is a support with sufficiently low resistance subject to strength,
coating and fabrication constraints. The numerical definition of the paper's
critical/matching line needs its simulation details before reproduction.
We should declare our own saturation criterion explicitly, for example a
specified remaining support share, rather than extract an exact optimum by eye.

Figures 5 and 6 make the intermediate layer important to the research design.
Low support resistance alone does not guarantee that the selective layer
becomes dominant: the intermediate layer may remain the main hydraulic load.
Report all three shares separately in a material-specific study.

## Traceable numerical information

All entries below are **author-model settings or reported model quantities**,
not measurements of a specimen available for calibration.

| Location | Quantity or setting | Interpretation |
| --- | --- | --- |
| Section 2.1 | Water viscosity 1.0e-3 Pa s | Assumed model viscosity; does not identify experimental temperature. |
| Fig. 2b caption | Particle diameter 3 µm; thickness varied over 1–9 mm | Single-layer support model and literature comparisons; exact experimental rows are in Table S2. |
| Fig. 4a–d caption and axes | Top-layer thickness 10 µm, particle diameter 1.7 µm; support particle diameter 2.8–3.2 µm and thickness 0.5–1.5 mm | Closest fixed-top-layer sweep. Panels e/f vary the stated top-layer setting. |
| Figs. 5–6 captions | Baseline bottom: 3 µm particles, 1 mm thick; intermediate: 1.7 µm particles, 10 µm thick; top: 0.3 µm particles, 1 µm thick | Three-layer reference configuration; not a reconstructed pore geometry. |
| Section 3.3, paragraph following Fig. 4 in the converted layout | Bottom resistance 1.01e10 m^-1; intermediate resistance 3.88e11 m^-1 | Values stated in prose. Their calculation must be reproduced before treating them as a quantitative model benchmark. |

**Particle diameter is not pore diameter or capillary radius.** The paper's
3 µm particles do not correspond to our synthetic 4/5 µm capillary radii.
Figure 4's horizontal diameter axis concerns particles. No digitised curve
points or error bars were extracted during this review.

Using only the two reported resistances gives a useful conditional bound.
For nonnegative top-layer resistance and unchanged intermediate/top layers,
complete removal of the bottom resistance can save at most

```text
Rs / (Rs + Rintermediate + Rtop)
    <= 1.01e10 / (1.01e10 + 3.88e11)
    = 0.0253705...  (about 2.54%)
```

This is an independently calculated **upper bound within that reported model
configuration**, at fixed flux, area and pump efficiency. It is not an
experimental energy-saving prediction. Any positive top-layer or other fixed
resistance lowers the bound further. The much larger savings in our earlier
synthetic sweep describe different assumed resistance partitions and remain
valid as conditional scenarios.

The corresponding L/R permeability conversions are approximately
9.90e-14 m² for the bottom layer and 2.58e-17 m² for the intermediate layer.
These are derived from the paper's reported model values, not measured OpenLB
permeabilities or independently established material properties. Production
parameters have not been replaced with them.

## Symbol and model boundaries

| Paper notation | Meaning | Our contract |
| --- | --- | --- |
| R1, R2, R3 | Bottom, intermediate and top hydraulic resistances, m^-1 | Support, separately identified fixed/intermediate, and selective resistance |
| alpha_i | Specific resistance R_i/L_i, m^-2 | Distinct from our dimensionless scenario ratio `alpha=Rselective/Rs,baseline` |
| k in Eq. 2 | Empirical modified Carman–Kozeny coefficient | Distinct from our Darcy permeability in m² |
| d_p | Constituent particle diameter | Cannot be passed as a capillary radius |
| r_m in Eq. 10 | Hydraulic average pore-size estimate from permeance, thickness and porosity | Does not directly identify selective-layer pore size or demonstrate unchanged rejection |

The paper's particle-packing/force model and its rejection model are additional
constitutive assumptions. Reproducing them would be a separate model task.
We should keep OpenLB as the detailed flow solver and use the paper as a
reference, rather than substitute its empirical permeability formula for the
required pore-scale stage. The present Julia layer model already accommodates
a supplied intermediate resistance; a numerical-model change is unnecessary
to preserve that distinction.

## Conversion checks before any equation implementation

The metadata reports OCR extraction on source pages 0–3 and text extraction
on pages 4–8. Several expressions in the Markdown need checking against the
original PDF. This is an extraction audit, not a claim that the published
paper has these errors:

- Eq. 1 labels the total denominator as R1 in one expression, while the
  accompanying prose and Eq. 4 use the sum of layer resistances.
- Eq. 2 displays a particle-size subscript differing from the subsequent
  definition. Its porosity exponent and empirical coefficient require checking;
  do not silently substitute a standard Carman–Kozeny formula.
- The denominator and sign rendered in Eq. 9 do not follow directly from
  differentiating Eq. 7 as supplied. Verify both original expressions and
  force-direction conventions before coding the balance.
- Eq. 18 renders a pressure symbol inconsistently. Preserve the independently
  defined permeance contract J/delta_p with its units.

Also check the meaning of the global hydraulic pore-size estimate before
using it in a rejection calculation. Holding top-layer structure fixed while
changing support resistance does not make a whole-membrane inferred pore size
an independent measurement of selective-layer pore geometry.

## Supplement received and next original source

The supplied `SupportingInformation.docx` contains S1, Tables S1–S3, Figure S1
and the complete supplementary reference list. All eight rendered pages and
all table cells have now been inspected. No Markdown conversion is needed
from the user. [The structured extraction](mo2022_supplement_extraction.json)
records six permeance and four pore-size entries in S2, with no row reporting
both. S3 contains nine wastewater comparisons. Table S1 gives model conditions,
not confirmed test conditions for those independent studies, and its density
unit needs clarification. No numerical matching-line tolerance or MATLAB code
is supplied. Main-equation typography still needs the original main PDF.

**Zhu et al. (2021), DOI 10.1016/j.ceramint.2021.01.088** has since been received
as full Markdown and figures. [Its matched-series study](../docs/zhu_matched_validation.md)
fits a limited M1 impedance model and tests held-out digitized points: nominal
maximum error is 4.53%, while conservative reading corners reach 10.78% against
a 10% project screen. Experimental uncertainty and connected geometry remain
unverified. No supplement is needed for this bounded result. [The supplement
review](mo2022_supplement_review.md) also records S2 source-scope ambiguities.

A collection of independent support measurements alone cannot identify
selective and intermediate resistances of another multilayer specimen.

The next geometry remains a graded support with fixed selective-side opening,
thickness and representative area, composed from stock OpenLB indicators.
Select its material scale and layer resistance partition after this evidence
check. No simulation source, production setting or completed result was changed
during the article review.
