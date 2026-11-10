# Literature and evidence

Use `sources/<paper-id>/` for a paper's metadata, single-paper extractions,
digitized points and review. The paper folders link shared records instead of
duplicating measurements. `reviews/` holds comparisons that involve several
papers; example-specific protocols and transferred scenario ledgers live with
their examples.

| Source | Role |
| --- | --- |
| [alhamed2009](sources/alhamed2009/README.md) | Audited date-stone carbon source; batch qualification and separate failed column comparisons. |
| [biesheuvel1999](sources/biesheuvel1999/README.md) | Support water-permeability/porosity/processing and bend-strength reference; distinguishes experimental Table 2 from model Table 3 and maximum-flow design curves. |
| [chiu2005](sources/chiu2005/README.md) | Experimental protocol and related ceramic/TiO2 series lead |
| [chiu2006](sources/chiu2006/README.md) | Concrete biological-water-treatment alternative with ceramic channels |
| [frederic2018](sources/frederic2018/README.md) | Support/selective-layer permeability contrast, spatial flux and fouling model |
| [garcia2008](sources/garcia2008/README.md) | Matched critical-flux/circulation-energy benchmark lead |
| [iheanacho2021](sources/iheanacho2021/README.md) | Alternative source candidate only; no extracted numerical dataset is used in the selected models. |
| [jan2025](sources/jan2025/README.md) | Whole-membrane experimental permeance and coating/interface context; not a separated-resistance calibration or fixed-selective-layer support comparison. |
| [meneghini2017](sources/meneghini2017/README.md) | Author dissertation; fitted coconut-shell parameters transferred to a synthetic scenario. |
| [meneghini2021](sources/meneghini2021/README.md) | Related publication; not independent validation of the dissertation parameters. |
| [mo2022](sources/mo2022/README.md) | Fixed-top-layer support plateau and explicit intermediate-layer resistance. Author-model values retained separately from experiment; no production parameter replacement. |
| [zhu2021](sources/zhu2021/README.md) | Matched support/transfer-coated whole-membrane M1 effective-impedance reference for bounded scale translation; no claim of exact specimen reconstruction. |
| [wakao_funazkri1978](sources/wakao_funazkri1978/README.md) | Existing packed-bed film-correlation reference |

## Shared reviews and original source ledger

- [Filtration evidence](reviews/filtration/evidence_review.md)
- [Crossflow paper comparison](reviews/crossflow/papers_review.md)
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
