# Internal transport at fixed catalyst mass

**Bottom line.** The largest bead change anchored to the available geometry
and the published reference effectiveness lowered conversion by at most
5.78 percentage points in the 2025 reactor calculation. The same case kept
all 20 predictions high, with 6.99 points RMSE. An extreme diffusion limit
brings the error down to 3.81 points and requires a much stronger reference
diffusion limitation. The operating-liquid transport contrast remains
unmeasured. [Run record](reference/summary.json)

## Evidence behind the scenarios

- The target study gives a 593 µm mean bead diameter. Its dry
  nitrogen-sorption estimate for particle porosity is 0.22 and its
  aqueous-tracer estimate is 0.40. Using their ratio gives an illustrative
  effective diffusivity ratio of 0.55 at fixed molecular diffusivity,
  tortuosity and constrictivity. The two porosity estimates come from
  different methods and fluids. The scenario assigns the smaller
  effective diffusivity to the packed bed. [Lorenzo et al. (2025)](https://doi.org/10.1021/acs.iecr.5c02780)
- Optical microscopy of commercial Novozym 435 on the same Lewatit resin
  used different bead samples from each sieve fraction. Equimolar ethyl
  butyrate and racemic 2-pentanol gave wetted mean diameters 9–21% above
  their dry group means. The largest sieved contrast, 256.2 versus 309.5 µm,
  gives a diameter ratio of 0.828 for the restricted-state sensitivity.
  The analogue also exhibited up to a factor-of-two initial rate change
  across the size classes. A wet-to-dry diameter difference in this liquid
  supplies a geometry comparison for the target catalyst. [Chaussard et al. (2024), SI Tables S1 and S2](https://doi.org/10.1021/acsengineeringau.4c00030.s001)
- Infrared mapping of Novozyme 435 found protein in an 80–100 µm exterior
  shell. The two shell calculations place the same amount of active
  catalyst in those layers at the reference diameter, with the inner
  radius fixed as a fraction of bead radius when the diameter changes.
  The mapping measured protein location. Activity and reactive-liquid
  wetting of the target beads remain unknown. [Mei et al. (2003)](https://doi.org/10.1021/bm025611t)

The combined porosity and diameter case deliberately stacks two
contrasts. Their simultaneous occurrence in this column is an
assumption. The [evidence record](evidence.json) retains the measurement
conditions and limitations.

A separate dry-porosity check uses the reported nitrogen pore volume
of 0.19 cm³/g and interprets the adopted polymer density of 1.17 g/cm³
as a skeletal density. Dividing the pore volume by pore plus solid
volume gives 0.182, and the corresponding diffusivity ratio is 0.455.
This calculation depends on treating the polymer density as the
skeletal density of the catalyst specimen.

## Calculation and held-out comparison

For a uniform spherical bead, the first-order modulus is
`Phi = R sqrt(k_vol / D_eff)` and the effectiveness is
`eta = 3 (Phi / tanh(Phi) - 1) / Phi^2`. Table B1 of the target study
calculates reference internal effectiveness values of 0.852–0.884 using
a borrowed tortuosity. Inverting the uniform-sphere expression gives
reference moduli of 1.44–1.68. These anchors represent the source model.
They leave the actual batch modulus uncertain.

The intrinsic rate per bead is fixed, so the intrinsic volumetric rate
scales with inverse bead volume. For this assumption,
`Phi^2` scales as `1 / (D_eff d)`. The shell solution sets zero flux at
its inner radius and preserves the amount of catalyst when its active
volume decreases. All scenarios use the ratio of restricted-state to
reference-state effectiveness at the same enzyme distribution. The ratio
scales both forward and reverse terms and leaves their equilibrium fixed.
The reactor still contains **6.0 g** for every condition.

The [2026 apparent batch kinetics](../../pentyl_propanoate/REPORT.md)
already reflect transport in the batch beads. They receive the relative
correction and retain the measured-water dispersion scenario. The
[2025 reconstruction](../../pentyl_propanoate_transfer/published/README.md)
uses the same relative correction around its Table B1 effectiveness and
retains its inlet convention and dispersion. Each model uses the same
reference coefficients and catalyst mass as its prior calculation.
Both forecast all 20 runs before opening the observations.

| Scenario | 2026 batch law RMSE (pp) | 2025 reconstructed law RMSE (pp) |
| --- | ---: | ---: |
| Reference | 11.84 | 10.46 |
| Porosity contrast only | 10.37 | 8.58 |
| Porosity and largest sieved diameter contrast | 9.77 | 7.83 |
| Skeletal-density dry porosity and largest sieved diameter contrast | 9.10 | 6.99 |
| Both contrasts, enzyme in the outer 80 µm | 10.68 | 8.99 |
| Both contrasts, enzyme in the outer 100 µm | 10.47 | 8.72 |
| Uniform sphere, limiting strong diffusion at the printed porosity contrast | 6.06 | 3.81 |
| Both contrasts and omission of all produced water in the sample | 8.81 | 6.87 |
| Skeletal-density porosity, diameter contrast and water omission | 8.16 | 6.08 |

The combined uniform-bead contrasts give apparent-rate multipliers of
0.867–0.890 with the printed dry porosity and 0.828–0.856 with the
skeletal-density check. Confining the same enzyme to the outer layers
reduces the penalty. In the strong-diffusion limit the
uniform-sphere multiplier tends to `sqrt(D_ratio * d_ratio) = 0.675`
for those same contrasts. This limit comes from `eta ~ 3 / Phi` as
`Phi` becomes large. `Phi * eta` is increasing for a uniform sphere,
so the limiting multiplier is a lower bound in that model. Its 2025
predictions have a 3.81-point RMSE, with 13 of 20 high. A reference
modulus of 1.44–1.68 is far from this limit.

For scale, attaining the earlier 0.64 exposure multiplier with that
reference range at fixed bead diameter would require a diffusivity
ratio of 0.132–0.160. This calculation uses the earlier volume
convention as its target. The held-out conversions score the fixed
scenarios.

At 60 °C and 3 mL/min, the 2025 result moves from 60.23% to 56.26%
under the two bead contrasts with the printed porosity. The skeletal-density
porosity gives 54.98%, and applying the water-omission readout as well
gives 53.44% versus 49.7% measured. The 2026 result with the same combined
sensitivity gives 63.30%.
[Conversion table](reference/bed_comparison.csv)

## Reproduce

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transport.particle.run
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_transport/particle/test_model.jl
```

The Julia controls compare the shell expression to independent radial
integration. The reacting runs compare two axial grids and conserve
normalized acid flux to 4.38 × 10⁻¹⁰. The largest grid change in
conversion is 0.00118 percentage points. The reference predictions
reproduce the preceding transfer audit.
