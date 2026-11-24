# Literature and evidence

Use `sources/<paper-id>/` for a paper's metadata, single-paper extractions,
digitized points and review. The paper folders link shared records instead of
duplicating measurements. `reviews/` holds comparisons that involve several
papers; example-specific protocols and transferred scenario ledgers live with
their examples.

| Source | Role |
| --- | --- |
| [alhamed2009](sources/alhamed2009/README.md) | Audited date-stone carbon source; batch qualification and separate failed column comparisons. |
| [bajic2017](sources/bajic2017/README.md) | Reactor assembly and nominal channel-volume convention; reported tubing diameter, no separate wetted fixture volume or narrow-channel quantitative pulse RTD. |
| [biesheuvel1999](sources/biesheuvel1999/README.md) | Support water-permeability/porosity/processing and bend-strength reference; distinguishes experimental Table 2 from model Table 3 and maximum-flow design curves. |
| [chiu2005](sources/chiu2005/README.md) | Experimental protocol and related ceramic/TiO2 series lead |
| [chiu2006](sources/chiu2006/README.md) | Concrete biological-water-treatment alternative with ceramic channels |
| [damstrup2007](sources/damstrup2007/README.md) | Independent continuous catalyst-inventory comparison, wet-layer expansion and matched-density scale-up; separate glycerolysis chemistry. |
| [frederic2018](sources/frederic2018/README.md) | Support/selective-layer permeability contrast, spatial flux and fouling model |
| [garcia2008](sources/garcia2008/README.md) | Matched critical-flux/circulation-energy benchmark lead |
| [iheanacho2021](sources/iheanacho2021/README.md) | Alternative source candidate only; no extracted numerical dataset is used in the selected models. |
| [jan2025](sources/jan2025/README.md) | Whole-membrane experimental permeance and coating/interface context; not a separated-resistance calibration or fixed-selective-layer support comparison. |
| [lerin2011](sources/lerin2011/README.md) | Direct reuse-cycle enzyme-activity assays distinct from conversion and recovered catalyst mass. |
| [lorenzo2025](sources/lorenzo2025/README.md) | Water-tracer and held-out pentyl propanoate conversion measurements, [batch markers](sources/lorenzo2025/batch_2to1_markers.md) and [published reactor model](sources/lorenzo2025/author_model.md). |
| [lorenzo2026](sources/lorenzo2026/README.md) | Complete 2:1 and 1:1 pentyl propanoate enzyme batch kinetics and separately reported apparent equilibrium; hot equimolar activity-loss challenge. |
| [vanhecke2024](sources/vanhecke2024/README.md) | Independent same-acid nonlinear inhibition shape and concentration scale, not a measured pentyl rate or deactivation constant. |
| [ascorbyl2019](sources/ascorbyl2019/README.md) | Qualitative Novozym 435 partial inactivation/retained activity evidence in different solvents. |
| [kuperkar2014](sources/kuperkar2014/README.md) | Same-acid isobutyl propionate loading curves, denominator unresolved. |
| [paroul2010](sources/paroul2010/README.md) | Geranyl propionate loading series and discordant nominally identical curve checks. |
| [richetti2010](sources/richetti2010/README.md) | High-loading palmitate early progress and caption rates. |
| [rejasse2004](sources/rejasse2004/README.md) | Independent time-resolved conventional-heating activity after substrate exposure; early exponential versus late held-out readings. |
| [sbardelotto2018](sources/sbardelotto2018/README.md) | Geranyl butanoate loading series and uncatalysed control. |
| [meneses2019](sources/meneses2019/README.md) | Fed-batch loading context; different acid exposure, not a fixed-charge rate test. |
| [meneghini2017](sources/meneghini2017/README.md) | Author dissertation; fitted coconut-shell parameters transferred to a synthetic scenario. |
| [meneghini2021](sources/meneghini2021/README.md) | Related publication; not independent validation of the dissertation parameters. |
| [mo2022](sources/mo2022/README.md) | Fixed-top-layer support plateau and explicit intermediate-layer resistance. Author-model values retained separately from experiment; no production parameter replacement. |
| [strnisa2018](sources/strnisa2018/README.md) | Six packed microreactor configurations with independent batch initial rate, water tracer and pressure measurements, and withheld geometry conversion series; synthetic packing supplement reviewed. |
| [panjan2018](sources/panjan2018/README.md) | Doctoral thesis Fig. 25 dimensional glucose pulses and separate detector characterization; likely shared blank/R1 observations, additional distinct middle-length pulse, incomplete wide-bed curves. |
| [zhu2021](sources/zhu2021/README.md) | Matched support/transfer-coated whole-membrane M1 effective-impedance reference for bounded scale translation; no claim of exact specimen reconstruction. |
| [wakao_funazkri1978](sources/wakao_funazkri1978/README.md) | Existing packed-bed film-correlation reference |

## Shared reviews and original source ledger

- [Filtration evidence](reviews/filtration/evidence_review.md)
- [Crossflow paper comparison](reviews/crossflow/papers_review.md)
- [Biocatalysis independent inventory and activity audit](reviews/biocatalysis_inventory/REPORT.md)
- [Original source ledger](sources/index.json)

The original JSON/CSV extractions and ledgers retain their exact bytes. Old
project-relative paths inside them are intentional provenance and resolve
through the [resource relocation registry](../docs/history/resource_organization.json).
Per-paper metadata provides the current navigation. Earlier access/status
statements remain dated evidence, not a replacement for a later audited record.

## Supplying articles

Place article Markdown, original PDFs, figures, conversion metadata and available
supplements in `inbox/<paper-id>/`. The inbox remains private and ignored by Git.
Keep source, specimen, figure/table, conditions, units and uncertainty with an
extraction. Unreported values stay unavailable or become explicit scenario
assumptions; measurements of unrelated specimens are not subtracted to infer
layer resistances. New experiments or inaccessible raw datasets are not gates.
