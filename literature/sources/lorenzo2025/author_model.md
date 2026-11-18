# Lorenzo et al. 2025 reactor model source audit

The [article](https://doi.org/10.1021/acs.iecr.5c02780) and its
[supporting information](https://doi.org/10.1021/acs.iecr.5c02780.s001)
give a steady, axially dispersed reactor model for pentyl propanoate production.
This record uses the 2:1 alcohol-to-acid feed. The PDF hashes and the Table S3
and Table B1 inputs are in [the model extraction](author_model.json). The article
says the reaction rate was based on batch experiments and reports the resulting
reactor calculations as solid lines in Figure 6c.

## Kinetics from the supplement

Supporting Information Equation S1 gives a Bi–Bi Ping-Pong reaction law for
immobilised lipase Lipozyme 435. It includes inhibition by both substrates, an
empirical acid denaturation term and a reverse term involving both products.
Write `A, B, W, E` for the acid, alcohol, water and ester concentrations.
Here `E` denotes ester. The numerator is
`C0cat * [k1(T) A B - k2(T) W E]`. The source defines the resulting `r`
as `-dC_acid/dθ`. The denominator is

```text
k3 A + k4 B + [k5 + Kd(T)] A B + k6 A² + k7 B² + k8 W
  + [k2_plain k3 k4 / (k1_plain k8)] E
  + [k9 + Kd(T)] A W + k10 B E
  + [k2(T) k4 k6 / (k1(T) k8)] A E + [k7 k8 / k4] B W
  + [k2(T) k4 [k9 + Kd(T)] / (k1(T) k8)
     - k2(T) [k5 + Kd(T)] / k1(T) + k8 k10 / k4] W E
```

Primes on the `k` coefficients are omitted here. The symbols
`k1_plain, k2_plain` reproduce the single-ester term's bare `k1, k2`.
The other terms show the temperature dependence explicitly. Whether the bare
symbols use reference or temperature-dependent values is unspecified.

```text
k_i(T) = k_i,ref exp[-k_i,exp (1/T - 1/T_ref)], i = 1, 2
Kd(T) = Kd_ref A0^Kd1 / (B0/A0)^Kd2 * exp[-Kd_exp (1/T - 1/T_ref)]
```

The propanoate parameter values are from Table S3. The table gives
`T_ref = 313.15` with absolute temperature defined in kelvins in the
nomenclature. The kinetic constants, the catalyst concentration
`C0cat` and the numerical time unit used to fit them have no printed unit
declarations. Equation S2 minimizes acid concentration residuals. The table
reports RMSE 0.228. Figure S3 plots concentration in mol/L and time in hours.
It specifies catalyst mass as 10% of initial acid mass. Figure S4 gives
different catalyst mass fractions at 40 °C. The printed captions interchange
the ratios for panels (b) and (c) in both figures. The in-panel label of (c)
and the concentration span identify the 2:1 propanoate series. The graphs
rise from zero and the right-hand scale gives percent conversion.

Reference 25 gives an earlier acetate model in mole fractions and a method for
estimating molecular diffusivity. The earlier acetate experiments use a 40 mL
batch and 10% catalyst mass relative to acid. The reference supplies no 2025
propanoate bed rate normalization or reaction-mixture properties.
[Lorenzo et al. 2023, Sections 2.3 and 4 and Appendix A](https://doi.org/10.3390/pr11061640)

## Packed-bed equations and input basis

The printed Section 4 equations are

```text
Eq 2   D_ax C″ - (u/ε_b) C′ + η_I r_obs = 0
Eq 3   ε_b u d_p / D_ax = 0.20 + 0.011 Re^0.48
Eq 4   C_acid(0) = C_0,acid; C_acid′(L) = 0
```

Derivatives are with respect to axial location. Equation 2 includes
axial dispersion, convection and internal mass-transfer effectiveness.
Equation S1 defines a positive disappearance rate for acid, whereas Equation 2
adds it to the acid balance. An acid-consuming model must use the negative of
the Equation S1 rate as the source of acid. A positive term is consistent with
a balance for reacted acid. Equation 4 explicitly fixes inlet concentration
and sets a zero outlet gradient. The authors describe these conditions as
Danckwerts conditions. That attribution leaves the printed inlet form as
a fixed concentration. The simulation uses a fourth-order Lobatto IIIA
boundary-value method according to the article.

The source specifies a 16 mm bore, 140 mm length and 6.0 g of catalyst.
The mean bead diameter is 593 µm. The interbead and particle porosities
used in the source calculations are 0.64 and 0.40 from aqueous tracers.
The 2:1 feed is 3.3 mol/L propanoic acid and 6.6 mol/L pentan-1-ol.
The reaction proceeds without added solvent and generates water. Reagents
were dried and the residual water concentration is unavailable. The
geometric bed volume is 28.15 mL. This gives an interbead liquid volume
of 18.02 mL with the stated porosity. Dividing this volume by the four
feed flows gives 25.74, 18.02, 9.01 and 6.01 min, which agree with the
rounded residence times in Section 2.3. This arithmetic supports use of
`u = Q / A_cross` and of `u/ε_b` as interstitial velocity. It leaves the
separate catalyst concentration in Equation S1 unresolved.

The same 6.0 g gives 213.15 g/L on a geometric volume basis or 333.05 g/L
on an interbead liquid basis. For a liquid-volume reaction balance, translating
a genuinely mass-specific catalyst rate uses 333.05 g/L. Inserting
213.15 g/L into that balance gives 64% of the catalyst exposure, the value
associated with 3.84 g in the liquid-volume calculation. The article and SI
provide no explicit mapping from the 6.0 g inventory to the `C0cat` factor of
Equation S1. The approximate 2 cm³ catalyst volume printed in Section 2.2
also differs from the 10.13 mL occupied fraction implied by the geometry
and aqueous voidage. The mass and volume statements retain their source basis.

Equation 3 is the Chung and Wen correlation as printed in this paper. Appendix A
separately reports water-tracer moments and uses its own boundary convention to
estimate `D_ax`. Section 4 directs use of Equation 3 in the reacting model.
The mixture density and viscosity at the reaction setpoints are unavailable.
Table S1 gives measured pure-component densities at 25 °C of 810.92 kg/m³
for pentan-1-ol and 988.07 kg/m³ for propanoic acid. These are pure
component properties at a separate temperature.

## Mass transfer treatment

Appendix B prints the following relations.

```text
B.1    η_E = (1 - Da_obs)²
       Da_obs = r_0,obs / (k_S a c_cat c_0)
B.2    Sh = k_S d_p / D = 1.17 Re^0.585 Sc^(1/3)
       Re = d_p u ρ / μ; Sc = μ D / ρ
B.3    η_I = r_0,obs / r_0 = [1/tanh(3φ) - 1/(3φ)] / φ
B.4    φ = (d_p/6) sqrt[r_0 / (D_eff c_0)]
       D_eff = D ε_p / τ
```

The printed Schmidt expression has dimensions of diffusivity squared if
`D` is a molecular diffusivity. The usual dimensionless definition is
`μ/(ρD)`. Table B1 gives `η_E = 0.999` at all five temperatures for the
2:1 propanoate feed. Section 4 treats external mass transfer as negligible.
Equation 2 includes `η_I`. The printed internal factors are 0.852, 0.861,
0.869, 0.878 and 0.884 at 40, 50, 60, 70 and 80 °C respectively. Appendix B
uses tortuosity 5.3 from a cited octyl-ester study and refers to the 2023
acetate paper for the diffusivity calculation. The propanoate liquid and
pellet diffusivities at individual setpoints and any axial variation of
effectiveness remain unspecified. The tabulated factors provide the reported
inputs for a replay.

## Figure 6c predictions

[The digitized CSV](author_conversion_2to1.csv) contains 20 model
predictions at the temperatures and flow settings listed in Supporting
Information Table S2, runs 65–84. The published markers in Figure 6c
are experimental, whereas the solid lines are the authors' calculations.
The experimental conversion observations used for checking predictions
are in [the existing measurement extraction](conversion_2to1.csv).

Figure 6c was rendered at 720 dpi from page 5 of the verified article PDF.
The crop and pixel calibration are in the model extraction. For each
temperature, colored stroke centers between flow settings were selected
in a small vertical window about the visually identified curve. The
first and last 24 px of each straight segment were omitted to avoid the
overlaid markers. Linear regression of stroke centers located the ends of
each segment. The two adjacent estimates were averaged at each interior
flow. Varying the color threshold and endpoint exclusion moved the result
by at most 0.06 percentage points. The largest mismatch between the two
adjacent segment estimates was 0.08 points after rounding upward.
The extractor retains a working half-range of 0.5 percentage points
for reading the graph, based on half the approximately 1-point stroke width.
The figure supplies no uncertainty for the model curves. Model conversions
are stored as fractions, rounded to three decimals. Coordinates are in
cropped image pixels. The sample IDs identify the corresponding experimental
settings.
