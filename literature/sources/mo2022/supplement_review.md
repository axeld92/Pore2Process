# Mo supplement evidence and the next validation dataset

Reviewed 8 November 2026. The supplied `SupportingInformation.docx` is the
supplement to Mo, Li and Yang (2022), DOI
[10.1016/j.watres.2022.118658](https://doi.org/10.1016/j.watres.2022.118658).
It contains Section S1, Figure S1, Tables S1–S3 and 19 references. The original
file is unchanged. All eight rendered pages, including the blank second page,
were inspected; the tables were also read directly from Word's XML.
[Structured extraction](supplement_extraction.json) retains original
table cells, source/render hashes, page locations, SI conversions and missing
values. Full text and a local Markdown reading extraction remain in the ignored
article folder. No original main-article PDF has been supplied yet.

The supplement resolves the short experimental citations. It does **not**
provide a matched support/selective-layer dataset by itself. The best next
paper is **Zhu et al. (2021), Design and optimization of ceramic membrane
structure: From the perspective of flux matching between support and membrane**,
DOI [10.1016/j.ceramint.2021.01.088](https://doi.org/10.1016/j.ceramint.2021.01.088).
Its publisher abstract describes similar transfer-coated membrane layers on
supports with different fluxes. That makes it the closest candidate for our
question; the full experimental series still needs inspection before acceptance.
[Publisher abstract and section snippets](https://www.sciencedirect.com/science/article/abs/pii/S0272884221001103).

## What S1 supplies and what remains unspecified

Table S1, rendered page 3, gives **author-model settings**: particle diameter
0.3–3 µm, layer thickness 1 µm–9 mm, pressure drop 1 bar, viscosity 1.0 mPa·s
and temperature 298.15 K. These do not establish conditions of the independent
experiments compiled in S2/S3. We have not assigned that temperature, viscosity
or pressure to those specimens.

The printed particle density is **3.7 kg/m³**, confirmed in both the XML and
rendered table. Its unit is suspicious for an alumina particle density.
It is retained verbatim and marked unusable without clarification; it has
not been silently converted to 3700 kg/m³ or used in a solver.

Figure S1, rendered page 4, shows the sequence from operating conditions,
particle size and layer thickness to force/packing calculations, layer
resistances, flux, inferred pore size and rejection. It gives a cycle count
symbol N and equations to evaluate, but supplies **no MATLAB source**, numerical
N, solver tolerances, root-selection rules or numerical saturation/matching
criterion. Table S1 does not resolve all empirical/force constants. The
converted main-article equations still require an original-PDF audit.
We therefore have not claimed reproduction of Mo's full model.

## Table S2 separates into two incomplete groups

All ten entries are transcribed as printed on rendered page 5. A tilde is
treated as unavailable, never as zero. Particle diameter, porosity, specimen
ID, temperature, measurement pressure, uncertainty and separated layer scope
are not tabulated. The generic material label is α-Al₂O₃; original specimen
compositions and geometry must be checked.

| Original reference | Reported thickness mm | Reported permeance LMH/bar | Reported average pore size µm |
| --- | --- | --- | --- |
| Qin 2020 | 6 | 3243 | unavailable |
| Qin 2017 | 4.7 | 7600 | unavailable |
| Guan 2017 | 2 | 15000 | unavailable |
| Zhu 2021 | 3 | 13140 | unavailable |
| Liu 2018b | 1.5 | 15800 | unavailable |
| Silva 2021 | 1 | 30000 | unavailable |
| Ha 2017 | 4 | unavailable | 1.1 |
| Gu 2021 | 5.9 | unavailable | 1.27 |
| Guesmi 2020 | 2 | unavailable | 1.35 |
| Qin 2015 | 3 | unavailable | 1.22 |

Thus, **six rows give permeance and four give pore size; none gives both**.
They are not ten fully specified specimens suitable for fitting a geometry
model. The pore-size values cannot be combined with permeances from different
rows. A model assumed particle size of 3 µm is not a measurement of every
compiled specimen. Permeance can be converted to SI, but deriving a verified
intrinsic k additionally needs the correct transport geometry and viscosity.

## Original-source checks already identify important distinctions

Publisher-indexed records and available methods/abstract snippets were checked.
Direct page fetches for Elsevier, ACS and RSC returned HTTP 403; these checks
do not amount to a full original-PDF audit. No institutional credentials were
used or requested.

- **Zhu 2021**, DOI above: the reported experimental design varies supports
  while retaining similar membrane-layer thickness and pore-size distributions.
  It also examines penetration resistance. Full data must establish specimen
  pairing, comparable top-layer properties, test conditions and uncertainty.
  Similarity in an abstract is not proof of identical selective resistance.
- **Qin 2017**, DOI
  [10.1016/j.ceramint.2016.09.206](https://doi.org/10.1016/j.ceramint.2016.09.206):
  the publisher's experimental snippet identifies bare tubular support
  permeance of 7562 LMH/bar, consistent with S2's rounded 7600. Its OD/ID are
  13.4/8.7 mm: radial wall thickness is 2.35 mm, while S2's 4.7 mm equals
  OD minus ID. The convention needs checking before using S2 thickness as
  a Darcy transport length. Coated-membrane permeance is a different quantity;
  no layer resistance has been inferred from the snippets.
  [Publisher](https://www.sciencedirect.com/science/article/pii/S0272884216317370).
- **Liu 2018b**, DOI [10.1039/C7RA13195J](https://doi.org/10.1039/C7RA13195J):
  original Table 1 identifies 15800 LMH/bar as a **bare support** property.
  Methods describe disks 25 mm across and 1.5 mm thick, whereas the table
  caption says tubular. The paper also reports coated membranes with much
  lower permeance. These specimen scopes cannot be interchanged. Its 25 °C
  gas-porometry condition does not establish water-permeance temperature.
  [Publisher methods and table](https://pubs.rsc.org/en/content/articlehtml/2018/ra/c7ra13195j).
- **Gu 2021**, DOI [10.1021/acsami.1c08157](https://doi.org/10.1021/acsami.1c08157):
  the abstract describes kaolin/alumina composite multichannel supports. S2's
  generic alumina label cannot identify which composition/specimen its row
  represents or establish a simple planar transport geometry.
  [Publisher](https://pubs.acs.org/doi/10.1021/acsami.1c08157).
- **Ha 2017**, DOI
  [10.1080/17436753.2017.1369658](https://doi.org/10.1080/17436753.2017.1369658):
  online in 2017, journal volume in 2018; compares coated alumina and
  diatomite–kaolin supports. The exact S2 specimen remains to be identified.
- **Guan 2017**, DOI [10.1111/ijac.12666](https://doi.org/10.1111/ijac.12666):
  publisher bibliographic identity is confirmed. It is a useful secondary
  priority for support microstructure, thickness and strength; its full
  experimental data have not yet been audited.

## Table S3 is a separate wastewater comparison

The nine rows on rendered page 6 include rejection, pore size, thickness,
permeance and references for different oils or organic dispersions. Their
figures/symbols were inspected and retained in the extraction. The table does
not give common temperature, pressure, concentration, droplet distribution,
crossflow, fouling state or uncertainty. Its permeances must not automatically
be treated as comparable clean-water values, and rejection is not a validation
of a hydraulic-only model. These are retained for a later separation showcase,
without changing our clean-water calibration or production settings.

## Concrete next step

**Zhu et al. (2021), DOI 10.1016/j.ceramint.2021.01.088**, has since been received
as full Markdown, metadata and figures in `literature/inbox/Zhu2021Design/`.
[The completed bounded study](../../../docs/verification/filtration/zhu_jpeg.md) digitizes
21 M1 points, fits odd ranks and tests even ranks using Julia. Nominal held-out
predictions pass a 10% project screen; the conservative reading-corner screen
fails. No supplement was identified or needed for this result, and another
broad article search is not needed. Mo's original main PDF remains
useful for equation typography, but the Word supplement is fully received and
does not need conversion or resubmission.

For Zhu, extract the support/paired coated-specimen series, test pressure and
temperature, effective area convention, normal-flow geometry, layer thickness
and pore distributions, and replicate/error information. First assess whether
the remainder after support resistance is sufficiently constant within a
matched series, including interface/penetration effects. Such a remainder is
an effective non-support resistance, not automatically a separately measured
intrinsic selective-layer property. Specify calibration/holdout groups and a
criterion before fitting, retaining uncertainty and checking pressure linearity
where data permit. Then propagate the measured resistance through the existing
Julia layers model and saved-data DWSIM boundary.

An OpenLB geometry comparison additionally needs specimen-linked structural
constraints and grid/domain/pressure validation. The previous Biesheuvel
held-out screen remains failed; it has not been replaced by accepting S2 rows
or retuning to this supplement. No simulation code, original OpenLB file,
dependency lock or completed result was changed in this source review.
