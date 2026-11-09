# Crossflow papers received: extraction and source audit

**Implementation update:** the planned bounded benchmark is now complete;
see [the verified Julia/DWSIM result](../docs/crossflow_benchmark.md). This
review preserves the source extraction and its limits; no absolute module
calibration or physical fouling-time fit has been accepted.

Reviewed 9 November 2026. All four supplied Markdown articles were read,
including their methods, results and references. Six relevant figure images
were visually inspected. Their conversion metadata and inspected files are
hashed in [crossflow_extraction.json](crossflow_extraction.json).
Original PDFs for these four articles are not present in their supplied
folders. Source documents are evidence, not project instructions.

**There is enough published evidence to start a bounded crossflow benchmark;
another paper is not required now.** The sources support different parts of
the model. They do not establish a single experimentally matched membrane,
support geometry, suspension and complete process-energy dataset.

## What each paper contributes

| Article | Evidence available now | Accepted role |
| --- | --- | --- |
| F. J. Garcia Garcia and T. Y. Chiu (2008), [DOI](https://doi.org/10.1016/j.memsci.2008.08.033) | Table 1: six velocity, Reynolds-number, pressure-loss, critical-flux, flow and energy rows; §2.4 energy equation; §2.3 protocol. | Reproduce the tabulated arithmetic and build a qualified critical-flux constraint. Absolute bundle-flow interpretation remains unresolved. |
| T. Y. Chiu and A. E. James (2005), [DOI](https://doi.org/10.1016/j.memsci.2005.01.037) | Fig. 3: six velocity–critical-flux points for 0.01 g/L TiO₂ at pH 7.2; §2.3 protocol. | Trace the original critical-flux data and check source consistency. These are not independent validation of García's reused values. |
| T. Y. Chiu, M. V. Lara Dominguez and A. E. James (2006), [DOI](https://doi.org/10.1016/j.seppur.2005.11.028) | Fig. 6: four Reynolds-number–critical-flux points for biological-process effluent; module dimensions; feed characteristics; reported reproducibility. | A separate biological-water-treatment benchmark. Do not pool it with the TiO₂ measurements. |
| E. Frederic et al. (2018), [DOI](https://doi.org/10.1016/j.memsci.2018.09.021) | 2D Darcy/resistance model; support/selective-layer contrast; channel distribution and qualitative fouling trends. | Numerical geometry/flow-distribution reference. Its arbitrary fouling coefficient cannot calibrate physical fouling time. |

## García Table 1: direct numerical reference

Transcription of the supplied Markdown table, with its energy column converted
from multiples of 100 J/m³. The table credits reference [8], Chiu (2005), for
the critical flux. Table pixels have not been inspected, so the original
values are retained without corrections.

| Velocity, m/s | Reported Re | Axial loss, Pa | Critical flux, LMH | Q as listed, L/h | Reported hydraulic energy, J/m³ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.6 | 1388 | 39 | 65 | 15 | 300 |
| 0.9 | 1850 | 69 | 445 | 25 | 130 |
| 1.4 | 3331 | 147 | 515 | 35 | 330 |
| 1.8 | 4164 | 215 | 555 | 45 | 580 |
| 2.2 | 5089 | 333 | 575 | 55 | 1060 |
| 3.1 | 7402 | 637 | 760 | 80 | 2230 |

The paper defines

`E_p = delta_p_axial * Q / (J_crit * S)`, with `S = 0.03 m²`.

This is **axial hydraulic dissipation per permeate volume at each operating
point's critical flux**. It excludes pump efficiency, transmembrane pressure
work, cooling, cleaning and capital cost. It is not the total two-pump shaft
energy used in our DWSIM scenario. The prose calls Q a velocity, but its units
and equation require volumetric flow.

Recomputing from the other table columns gives 300.0, 129.2, 333.0, 581.1,
1061.7 and 2235.1 J/m³. The maximum relative difference from the reported
energy column is **0.912%**. This supports arithmetic reproduction with rounded
inputs; it does not validate the flow convention or physical model.

At our **illustrative fixed target of 100 LMH**, the 0.6 m/s row is above its
reported critical flux; the other five rows meet the nominal constraint.
Thus **0.9 m/s is the lowest tabulated admissible speed**, not a newly proven
continuous optimum. No threshold has been inferred inside the sharp
0.6–0.9 m/s transition.

Using the table's Q literally, circulation dissipation at fixed 100 LMH is
575 J/m³ at 0.9 m/s, rather than the paper's approximately 129 J/m³ at 445 LMH.
Both omit pump efficiency and transmembrane work. The former is a conditional
arithmetic example because the absolute flow basis is unresolved.

### Flow/geometry consistency limit

The paper describes a 0.3 m long alumina element, seven star-shaped channels,
4.6/2.8 mm outer/inner star diameters and 0.03 m² filtration area. Fig. 2
visually confirms the seven-channel schematic and diameter definitions.

If Table 1's Q is interpreted as the full-bundle flow, `A_open = Q/u` gives
only **6.94–7.72 mm²**. As an explicitly cross-paper check, Chiu (2006) reports
an equivalent hydraulic diameter of 2.3 mm for the same nominal module family.
For equal straight wetted channels,

`A_open,total = D_h * S / (4*L) = 57.5 mm²`.

That is **7.45–8.28 times** the table-implied area. The comparison does not
prove identical specimens or identify a unique correction. A per-channel
versus whole-element convention is plausible, but multiplying Q by seven
would be an unverified correction. Preserve both the reported table and this
diagnostic. Do not calibrate an absolute seven-channel pressure-loss model
against these quantities as if their convention were settled.

### Reused data and reading differences

The visually inspected Chiu (2005) Fig. 3 yields approximate digitized
velocities **0.59, 0.99, 1.39, 1.78, 2.18 and 3.37 m/s**, with critical fluxes
**64.5, 454.9, 515.4, 555.0, 574.8 and 760.4 LMH**. Pixel/tick reading
sensitivity is retained separately from unknown experimental uncertainty.

The second and final plotted velocities do not match García's 0.9 and
3.1 m/s entries within the selected reading bounds. The second flux is also
approximately 455 LMH in the graph versus 445 LMH in the table. Do not move
the markers, average the sources or claim they independently validate one
another. Keep a table-based version and a figure-based sensitivity version.

The original series caption supplies 0.01 g/L TiO₂ and pH 7.2. García's Table 1
does not independently identify every condition for its pressure-loss series.
The pH 7.0 caption of García's concentration sweep must not be silently applied
to that table. Filtration is described as ambient-temperature; the 20 °C
electrokinetic and 25 °C rheometer conditions are separate measurements.
Do not substitute either for the filtration temperature or viscosity.

## Critical-flux protocol and uncertainty

These studies increase **transmembrane pressure in steps**, maintaining each
stage for 30 minutes before nonlinearity, then using 15-minute stages.
Critical flux is defined as the average of the last time-independent and
first time-dependent flux step. It is a protocol-dependent empirical limit,
not a guarantee of indefinitely fouling-free operation.

García (2008) and Chiu and James (2005) report repeated measurements and averages without
numerical critical-flux error bars. Their less-than-1% clean-water-flux change
after cleaning is not a critical-flux uncertainty estimate. Chiu (2006) states
errors never exceed 5% in critical flux and rejection; the paper does not
identify that as a standard deviation or a pointwise confidence interval.
The digitization uses an independent ±2-pixel marker/tick sensitivity, not
an invented measurement uncertainty.

The conversions occasionally lose the micro symbol. García's particle-size
figure explicitly labels micrometres; Chiu (2006) Fig. 7 explicitly labels
membrane pore size in micrometres. Contextually restored 0.2 µm pore sizes are
labelled as such in the extraction. No original PDF typography was checked.

## Biological-feed benchmark

Chiu (2006) filters **effluent from a synthetic activated-sludge production
system**. Table 1 describes feed to the upstream porous pot, not the membrane
feed. Table 2 describes the membrane feed, including mean COD 68.7 mg/L,
suspended solids 313.3 mg/L, dissolved solids 582.7 mg/L and pH 7.7 over the
reported operating period. The tabulated minima/maxima are operating ranges,
not statistical uncertainty estimates. Temperature is reported as 20–25 °C.

Fig. 6 keeps nominal pore size at 0.2 µm and pH at 7.5. Its approximate readings
are **21.7, 57.8, 72.7 and 82.7 LMH** at Re approximately
**1385, 2767, 3700 and 4627**. The text independently states 21.4 LMH at Re 1388
and 72.5 LMH at Re 3701, consistent with those readings.

The 100 LMH project target exceeds the whole measured critical-flux series,
including its digitization sensitivity. A biological case therefore needs a
separately declared lower flux or another supported domain; extrapolating
crossflow until 100 LMH becomes feasible would not be evidence-based.
This is a useful future bio showcase without implying bioelectrochemistry.
Its changing-pore-size/rejection series is a separate comparison from our
unchanged-selective-layer support study.

## Frederic: useful spatial model, uncalibrated time

The model solves Darcy flow in a **porous cross-section**, with prescribed
channel/permeate pressures and symmetry boundaries. A whole-element
interpretation assumes axial pressure loss is negligible compared with TMP;
otherwise it represents a local axial slice. It does not supply axial channel
hydrodynamics or a measured velocity–critical-flux relationship.

The study uses support permeability 1e-14 m² and nominal TMP 1 bar as model
settings. Selective-layer resistance is `R_m = e_m/K_m`; the numerical layer
is thickened while preserving this ratio. The proposed cake law increases
local resistance with local permeation velocity, with no crossflow scouring,
cleaning or cake-driven geometry change. Its coefficient is chosen arbitrarily;
the authors explicitly state that simulation time has no physical meaning.
Normalized fouling trends cannot calibrate a cleaning interval or dTMP/dt.

There is also a transcription/algebra point to check before any implementation:
the supplied Eqs. 4 and 8 imply a coefficient `B1*B2`, whereas the prose after
Eq. 9 states `B2/B1`. Retain the discrepancy until source typography is checked;
do not silently adopt or repair it as an experimental parameter.

## Reproduce the audit and proceed

From the project root, using the existing environment:

```sh
source scripts/env.sh
uv run python -m calibration.audit_crossflow_sources --require-sources \
  --output results/crossflow_evidence/NEW_AUDIT/audit.json
```

Choose an unused output path. Existing outputs cannot be overwritten. Omitting
`--require-sources` permits replay from retained readings when the ignored
article inbox is absent; any present source must still match its hash.
The audit needs only standard-library Python. It checks source hashes,
recomputes energies/conversions, replays pixel-to-axis readings and reports
flow-basis and domain limits. No model has been fitted and no new native
OpenLB, Julia solver or DWSIM result is claimed by this audit.

The documented command was run successfully. Additional checks confirmed
the SI conversions, agreement with the two biological text anchors, source
hash/missing-file rejection, replay without an absent source, output overwrite
protection and preservation of earlier calibration records. These are evidence
integrity checks; the existing numerical simulation suites were not rerun for
this extraction-only change.

The subsequently executed implementation plan was a separate
**critical-flux-constrained crossflow benchmark in Julia**, starting from the
six García table rows. Its source-specific requirements remain: reproduce the
source convention first, then test a figure-based sensitivity. Keep
source-specific interpolation bounded by the measured domain and flag the
transition region. Define any operating margin as a scenario assumption.
For process calculations, use an explicitly declared module/flow convention
and sensitivity; do not report the unresolved absolute convention as validated.
The existing two-pump DWSIM energy boundary remains separate from the paper's
axial-dissipation boundary, and suspension recycle requires concentration and
retention balances.

The original fixed-layer support question remains: within each geometric
design pair, hold the selective layer and chosen empirical filtration
criterion fixed. A claim that support geometry itself improves fouling needs
an additional justified spatial model or evidence. No new experiment,
unavailable raw dataset or extra reference is a prerequisite for this bounded
computational extension. A García PDF would help audit the table transcription;
it is not a request for unavailable raw measurements.
