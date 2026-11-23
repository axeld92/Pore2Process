# Test protocol: catalyst mass proportionality

Use the frozen [published 2:1 BSTR measurements](../../../literature/sources/lorenzo2026/README.md)
and the four loadings at 40 °C, assays 1–4. Do not estimate
any rate or loading parameter from the 20 continuous-reactor
conversions. [Machine-readable design](protocol.json)

## Experimental information

Lorenzo and Ortega (2026), main Table 1 and SI Table S2, specify
2.5, 5, 10 and 20% catalyst relative to initial acid mass and
nominal 29.7 mL pentanol plus 10.3 mL propanoic acid in the
same 40 mL charge. Compute the approximate absolute catalyst
mass from the nominal acid volume and measured pure-acid
density at 25 °C, and the concentration per nominal initial
liquid volume. No operating-temperature slurry volume or
measured solids holdup is available.

The paper reports duplicate runs and a conversion uncertainty
of ±0.002, but does not supply the individual replicate curves.
The printed conversion and composition-derived conversion
disagree at some times. Retain printed conversion as the
primary response, and separately use
`1-x_acid/x_acid,0` and `x_ester/x_acid,0` as sensitivity
responses. Use the nominal 2:1 ratio by default and compare
the printed initial composition ratio as a separate
sensitivity. Section 2 reports liquid sampling at 0.1 mL
per sample and catalyst reuse in multiple assays; neither
assay-specific reuse history nor catalyst loss during
sampling is specified. In a liquid-only sampling screen,
reduce liquid inventory by 0.1 mL after each positive-time
measurement, retain all catalyst, and integrate exposure
with the resulting volume fraction. This bounds the scale
of that particular assumption without adjusting the data.
[Article](https://doi.org/10.1021/acs.iecr.6c02375)

## Direct test independent of a kinetic equation

At fixed temperature and feed, if the total rate is proportional
to catalyst mass, changing catalyst mass at fixed charge only
changes the reaction clock:

```text
exposure = (m_cat / n_acid,0) * time  [kg s/mol]
scaled time = time * (loading_pct / 10)
```

Plot all four unadjusted conversions on both ordinary and
scaled time. Use linear interpolation between measured
points to obtain the first time at 10, 20, 30 and 40% acid
conversion, levels reached by all four runs. Do not extrapolate
to higher conversions. For the primary response, retain the
bracketing measurement times and the crossing times at target
±0.002, a reading sensitivity rather than an uncertainty
distribution. The average rate per gram up to the same
conversion is inversely proportional to the exposure; the
average specific rate within a conversion interval is
`(X_b-X_a)/(exposure_b-exposure_a)`.

As a descriptive check, estimate a total-rate exponent
`alpha` from the slope of -log(time to chosen conversion)
versus log(loading), with independent intercepts per conversion
level. It is fitted using the specified training curves alone.
For each held-out loading, align the other curves by time
times (loading/10)^alpha, interpolate their conversions at
the held-out normalized times, and average the three training
curves. Compare alpha=1 and the training-derived exponent
only on the held-out points inside *every* training curve's
exposure range for both laws. This model-free assessment cannot
extrapolate to late equilibrium data.

## Reversible kinetic comparison

Use the earlier Julia apparent reversible second-order law
and its dry-feed/mole-fraction closure. Generalize to R, the
initial alcohol-to-acid molar ratio:

```text
dX/dt = (m_cat/n_acid,0) * k_ref * f(C/C_ref) *
        ((1-X)*(R-X) - X^2/Kx) / (1+R)^2
```

Here `k_ref` has units mol/(kg s), `Kx` is dimensionless
and incorporates the effective liquid nonideality, and
`f` is a relative apparent rate per gram at fixed conversion.
For the primary fit at R=2 the exact solution reproduces the
existing Julia kinetics model. The generalized solution uses
the gap `g=sqrt((R-1)^2+4*R/Kx)` and
`s=1-exp(-g*k_ref*f*exposure/(1+R)^2)`, yielding
`X=2*R*s/(2*g+(1+R-g)*s)`.

Compare three shared, separable rate laws, each normalized
to the 10%-loading catalyst concentration:

| Law | Relative rate per gram, u = C/C_ref | Interpretation |
| --- | --- | --- |
| Proportional | 1 | Total rate linear through zero |
| Power | u^(alpha-1), 0 <= alpha <= 3 | Sub- or superproportional total rate |
| Crowding | 1/(1+s*(u-1)), 0 <= s <= 1 | Nonincreasing rate per gram; total rate saturates |

The crowding law is a conditional engineering functional form,
not a derivation from measured bead interactions. At s=0 it
is exactly proportional; at s=1 total rate is constant for
positive loading. Both new models introduce just one shape
parameter. For each shape, fit k_ref and a common Kx using
the training conversions. Profile the bounded shape using
a uniform scan (25 points) followed by deterministic
bracket refinement and explicit comparison to the nested
proportional model. Unresolved extreme-shape nuisance trials
can fall back to a bounded log-rate/log-equilibrium profile.
Check selected solutions converge in the interior and
repeat the full primary shape scan at 49 points.

## Holdout boundary and sensitivity

Hold out one complete loading (all its points), train with the
other three and repeat for each loading. Fit common equilibrium
and kinetic parameters anew in every fold: main Table 2
provides Kx=4.53 from the 10% experiment, so fixing that value
in the strict 10%-holdout fold would leak plateau information.
Use the fixed published Kx only for all-data and individual
loading diagnostics. The latter are descriptive fits, never
used to predict an individual held-out curve.

Equal weights per observation form the primary fit. Compare
equal total weight per batch series, printed initial ratio,
and the two reconstructed conversion responses separately.
Report the RMSE per held-out loading and its equal-loading
aggregate (square root of the mean of fold MSEs). The
predeclared practical screen asks for at least a 10% reduction
in equal-loading RMSE and no fold worsening by more than
0.2 percentage points. Interpret that threshold as resolution
context, not statistical significance. There are just four
loading conditions with unavailable replicate curves.

Test the closest measured pair (10% and 20%) separately,
training only on these and scoring the omitted 2.5 and 5%
curves with no refitting. Their two-point concentration
trend cannot distinguish the two possible one-parameter
loading-law forms outside those concentrations. Retain both
extrapolations as explicitly unqualified diagnostics even
when the lower-loading predictions deteriorate.

## Transfer beyond the measured concentration range

Compute the nominal original batch concentrations in g/L
and bed mass per interbead liquid volume. Use bed mass per
total cartridge volume as a separate proxy for a putative
concentration-sensitive interaction. This changes only the
argument of f, never catalyst inventory or reaction balance.
Use the 40 °C full-fit and high-pair rate factors at both
concentration proxies; extrapolate that same relative factor
across 50–80 °C only as a sensitivity because temperature
interaction with loading is unmeasured.

1. The 2026 reversible model uses its independently fitted
   10%-only rate and published equilibrium for each temperature,
   the measured-water nominal rectangular-pulse axial-spread
   sensitivity (Pe from the previous transfer audit), and 6 g
   catalyst. Apply the relative loading law only to k.
2. The 2025 published-rate reconstruction retains the preceding
   Wilson–Geankoplis film, internal effectiveness, geometry,
   inlet/outlet and dispersion conventions. In the local
   film/source calculation apply the same independent relative
   mass-specific rate multiplier; the film balance, actual
   catalyst concentration and solid/liquid inventories remain
   unchanged. This is a continuity sensitivity with a different
   rate law, not an additional fit.

Verify direct Julia kinetics against the earlier model and an
independent numerical integration at a non-2:1 feed, recover
known synthetic law parameters on withheld loading data, and
verify exact proportional replay and film balance. Compare
axial grids and conservation on every reactor sensitivity.
Write the frozen predicted results before opening the 20
experimental observations. Compute scores for all named
scenarios without using any score to choose a loading law.
