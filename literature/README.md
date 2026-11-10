# Article workspace and experimental evidence

## Qualified operating-map audit — 9–10 November 2026

[The operating-window ledger](operating_window_evidence.json) now audits the
four supplied crossflow papers for the higher-Re hydraulic map. Chiu (2006)
is the biological-feed critical-flux source; Chiu (2005)/García (2008) retain
their separate TiO2 scope and shared-data caveat, and Frederic (2018) supplies
no calibrated physical fouling time. No additional article was required for
this bounded milestone.

[The computational report](../docs/operating_window.md) distinguishes the
source's short-step protocol, figure-reading sensitivities, unavailable
pointwise uncertainty and transferred retention from verified native/process
hydraulics. Linear interpolation, monotonic bounds, equal-Re/equal-velocity
mapping and transfer factors remain separate scenarios. An internal
leave-one-point-out shape check is not independent experimental validation.
The synthetic geometry, Zhu impedance and biological feed do not constitute
a matched experimental specimen.

## Crossflow performance search: 9 November 2026

[The initial crossflow search](crossflow_evidence_review.md) identified eight
further leads after screening the original four main articles and Mo's
supplement. **García 2008, Chiu 2005, Chiu 2006 and Frederic 2018 have now all
been received** as full Markdown, metadata and figures and reviewed.
[The received-paper report](crossflow_papers_review.md) and
[extraction](crossflow_extraction.json) retain six García table rows and ten
Chiu figure readings. [The screening ledger](crossflow_screening.json) records
the initial access state and subsequent receipt, source hashes and limits.

The source audit reproduces the reported energy arithmetic within 0.912%.
It also identifies reused measurements and an unresolved absolute bundle-flow
convention. No new model has been fitted or accepted as physical calibration.
The [bounded critical-flux-constrained benchmark](../docs/crossflow_benchmark.md)
now passes 846 Julia assertions, 50 Python tests and 32 actual DWSIM scenarios;
absolute module calibration remains unresolved. Another article is not required
for that completed milestone. The biological-feed series is retained
as a separate option and does not support the existing 100 LMH target within
its measured critical-flux range. The four original PDFs are absent; a García
PDF could help audit table typography, without making raw data or new
measurements a prerequisite.

## Biological-feed sensitivity in the explicit channel

[The additional Chiu 2006 extraction](chiu2006_bio_extraction.json) retains
the membrane-feed COD of 68.7 mg/L (Table 2; reported operating range
65.3–72.1 mg/L), below-critical COD rejection of 0.65 (Table 3; 0.2 µm,
1.6 m/s, 0.1 bar, pH 7.5), and the lowest critical-flux point of 21.4 LMH
at Re=1388 (Section 3.2/Figure 6). It records conditions, source hash,
unavailable point-specific uncertainty and the Markdown unit-typography limit.

The [explicit-channel showcase](../docs/star_channel_showcase.md) uses a
different, independently verified low-Re hydraulic range. Its 0.5–0.8 constant
local COD-rejection range is a transferred sensitivity, not a fitted or
validated law at the new speed. The source critical-flux point remains context;
critical/sustainable flux and fouling time for the new channel remain
unavailable. COD is balanced externally in Julia; DWSIM uses water-carrier
properties. No new entry is added to the accepted physical-calibration set.

## Supplied material and completed calibration

Place article Markdown, original PDFs and supplementary data in `inbox/`.
The folder contents are ignored by Git so institutionally supplied full text
is kept locally. All three requested main articles, their available figures and
conversion metadata have now been received in `inbox/Mo2022Dissecting/`,
`inbox/Biesheuvel1999Design/` and `inbox/Jan2025Silicon/`. Use a subfolder per paper,
for example `inbox/mo_2022_multilayer/`, with `article.md`, `article.pdf` and
`supplement/`. Markdown should retain figure/table captions, equations, units,
references and page or section markers. Keep original plots when numerical
values only appear graphically. No credentials belong in this folder.

The selected application is **clean-water filtration through an asymmetric
ceramic membrane**. Alumina is the first literature candidate; SiC is retained
as a comparison. A limited alumina support reduced model has now been fitted
to Biesheuvel data and tested against held-out cases; **the 25% screen fails**.
See [the calibration report](../docs/experimental_calibration.md) and
[commands](../calibration/README.md). The separate OpenLB study still uses
synthetic geometry and assumed layer resistance. Exact reconstruction of a
manufactured specimen is outside [the computational scope](../docs/verification_scope.md).

The later [fixed-opening OpenLB comparison](../docs/graded_support_study.md)
uses Zhu's effective non-support impedance and a common support baseline
calibration. It now connects three measured geometric resistance ratios to
nine actual DWSIM pump cases. This is a conditional literature-informed
scenario; no additional experimental evidence was required to execute it.

[The combined evidence review](filtration_evidence_review.md) assigns a distinct
role to each paper. Biesheuvel Table 2 provides nine experimental
support-property cases, now accepted for this limited calibration; Jan provides a whole-membrane comparison; Mo
provides the fixed-top-layer model. These are not a matched three-paper
resistance dataset. Structured values and source hashes are retained in the
linked extraction records.

## Mo et al. main article and Word supplement received

**J. Mo, X. Li and Z. Yang (2022), “Dissecting the structure-property relationship
of ceramic membrane with asymmetric multilayer structures for maximizing
permselectivity,” Water Research 220, 118658.**
[DOI](https://doi.org/10.1016/j.watres.2022.118658),
[verified abstract and metadata](https://pubmed.ncbi.nlm.nih.gov/35640511/).

The supplied main article, all seven scientific figures and the Word
supplement have now been reviewed. [The review](mo2022_review.md) maps the paper to our model, identifies
the intermediate-layer contribution and distinguishes model quantities from
experimental evidence. [The extraction record](mo2022_extraction.json) retains
source hashes, units and locations. No Mo experimental dataset has been accepted.

The supplement is in `inbox/Mo2022Dissecting/SupportingInformation.docx`;
there is no need to resubmit or convert it. [The supplement report](mo2022_supplement_review.md)
and [extraction](mo2022_supplement_extraction.json) retain all ten S2 entries,
nine S3 entries, model settings and identified references. S2 has no row
containing both permeance and pore size; original specimens/conditions still
are recorded as limits on using those rows. Equation typography can be revisited
if that Mo model is implemented; this is not a prerequisite for the hydraulic case.

The subsequent article **W. Zhu et al. (2021), Design and optimization of
ceramic membrane structure: From the perspective of flux matching between
support and membrane**, [DOI 10.1016/j.ceramint.2021.01.088](https://doi.org/10.1016/j.ceramint.2021.01.088).
has now been received under `inbox/Zhu2021Design/` as full Markdown, metadata
and figures, followed by the original PDF. [The PDF refinement](../docs/zhu_pdf_refinement.md)
digitizes the same 21 M1 points from its larger native image and calls Julia
with the original model, split and 10% screen: nominal held-out maximum error
4.56%, conservative reading-corner maximum 7.88%. The earlier JPEG reading
(4.53% nominal, 10.78% corner) is retained. No supplement was identified or
needed. Unreported conditions/uncertainty stay unknown; no new dataset is
required before the synthetic geometry extension.

## Biesheuvel and Verweij: limited support calibration completed

P. M. Biesheuvel and H. Verweij (1999), “Design of ceramic
membrane supports: permeability, tensile strength and stress,” Journal of
Membrane Science 156, 141–152.
[DOI](https://doi.org/10.1016/S0376-7388(98)00335-4).
The supplied main article, metadata and all nine scientific figures have been
reviewed. Table 2 was transcribed with its original units and SI conversions
into [support_membrane_extraction.json](support_membrane_extraction.json).
It contains measured water permeability, porosity and bending-derived strength
for nine processing conditions. The same study's Table 2, formula and powder
sizes have now been visually checked in Chapter 7 of the public
[author thesis](https://ris.utwente.nl/ws/portalfiles/portal/6077974/t0000015.pdf).
[The audit](biesheuvel_source_audit.json) retains source/page/image hashes;
this is the same dataset, not another independent experiment. The journal
publisher PDF was not inspected. The one-parameter support closure now has a
reproducible fit and held-out check; its primary screen fails. Existing native
simulation parameters have not been replaced. Water-test temperature, raw
series, uncertainty and specimen-linked connected morphology are unreported
limits, not acquisition requirements.

## Evidence already found

The public full text of A. Jan et al. (2025), “Silicon carbide coated alumina
tight-ultrafiltration membrane prepared by low-pressure chemical vapor
deposition for sulphate ion retention,” Desalination 613, 119085 is available:
[DOI](https://doi.org/10.1016/j.desal.2025.119085),
[university-hosted PDF](https://ris.utwente.nl/ws/portalfiles/portal/501138486/1-s2.0-S0011916425005612-main.pdf).
It is a potential whole-membrane reference, not a separated-resistance dataset.
The supplied local Markdown and available figures have now been reviewed.
The uncoated membrane already contains a selective layer, and SiC deposition
changes that layer. Unreported coating/interface and hydraulic conditions limit
what can be inferred; they do not block our selected Zhu-based comparison. Details and
caveats are in [the combined review](filtration_evidence_review.md) and
[the source ledger](sources.json).

## Extract when reported

Use the available papers and supporting information only. The fields below
guide extraction; they are not a demand for a complete raw dataset or new
experiments. Store unreported fields as unavailable, and label any scenario
assumptions separately. For each relevant specimen, record what is supplied:

- Material, specimen ID, support thickness, selective/intermediate-layer
  thicknesses, pore-size distributions, porosity, dimensions and normal-flow direction.
- Pure-water flux versus transmembrane pressure, water temperature and how
  permeance was temperature-corrected, wetting/conditioning and test mode.
- Bare-support permeance and coated-membrane permeance for the **same support
  design and comparable conditions**, with uncertainties and replicates.
- Evidence that selective-layer microstructure, defect rate and rejection
  remain comparable when support geometry changes; interface infiltration matters.
- Mechanical strength or pressure rating, if comparing feasible material designs.

An extracted datum must retain DOI, page/section/figure/table, specimen scope,
value/unit, SI conversion and whether it was tabulated, digitised or inferred.
Keep reported conditions/uncertainty and explicit unknowns. Record extraction
method and distinguish reading sensitivity from measurement uncertainty.
Preserve published values before any temperature correction. A fitted
datum cannot also count as independent validation.

For the same fully wetted specimen under matched conditions, an apparent extra
resistance may be obtained from coated minus bare-support resistance. That
difference includes selective, intermediate and interface effects unless the
experiment identifies them separately. Never subtract values from unrelated
papers/specimens and call the result a measured selective-layer resistance.

The ledger accepts Biesheuvel's support table and Zhu's matched M1 figure for
their respective limited calibrations. Neither is a specimen-resolved 3D
geometry dataset. Failed model screens remain recorded alongside successful
software execution and other numerical checks.
