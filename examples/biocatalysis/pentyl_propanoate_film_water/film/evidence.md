# Liquid-film transport inputs

The screening uses pentan-1-ol as the transport medium and propanoic acid as the transferred solute. Pentan-1-ol makes up two thirds of the initial feed on a molar basis. Transport measurements for the reacting acid, alcohol, ester and water mixture at the five reactor temperatures were unavailable. The property scenarios use independent measurements and estimates. The reactor conversions were kept held out. [Lorenzo et al. (2025), reactor conditions](https://doi.org/10.1021/acs.iecr.5c02780)

## Density, dynamic viscosity and bulk diffusivity

Density and viscosity come from measurements of pure pentan-1-ol at 0.100 MPa. The 40, 60 and 80 °C values are direct table entries. The two intermediate values interpolate density and log viscosity linearly in temperature. Moreau et al. report expanded measurement uncertainties of 0.70 kg/m³ for density and 1.5% for viscosity, with coverage factor 2. Replacing the reacting mixture with neat pentan-1-ol introduces a separate, unquantified composition effect. [Moreau et al. (2021), Tables 4 and 9](https://uvadoc.uva.es/server/api/core/bitstreams/1dbed714-67ce-4929-a858-2ef042e4aba1/content)

The central estimate uses Wilke–Chang for dilute propanoic acid in pentan-1-ol with the measured viscosity and a Le Bas boiling-point molar volume of 90.6 cm³/mol for the acid. The solvent association factor of 1.0 is an assumption. The high scenario uses 1.5 so that the calculated diffusivity rises by √1.5. These estimates omit concentration dependence at the 3.3 mol/L acid feed. [Wilke and Chang (1955)](https://doi.org/10.1002/aic.690010222), [La-Scalea, Menezes and Ferreira (2005), Appendix Table A1](https://doi.org/10.1016/j.theochem.2005.05.030)

The low scenario uses neat pentan-1-ol self-diffusion as a slower analogue. Harris et al. measured it by NMR over −20 to 91 °C. Their fit has a 0.9% relative standard uncertainty for that liquid. The difference between neat-solvent self-diffusion and acid transport in the feed remains unknown. The low and high inputs are sensitivity scenarios with no assigned probability. [Harris, Ganbold and Price (2015), Equation 9](https://doi.org/10.1021/acs.jced.5b00246), [NIST ThermoML data](https://trc.nist.gov/ThermoML/10.1021/acs.jced.5b00246.html)

Two homologous 1-butanol systems provide literature context. At 60 °C a 1:1 acid-to-1-butanol feed had a dynamic viscosity of 0.931 mPa s and density of 869.9 kg/m³. Liu et al. measured propanoic acid diffusivity in water-saturated 1-butanol at 25 °C. The numerical table was unavailable in the accessible copy. [Lilja et al. (2005), Table 2](https://doi.org/10.1016/j.cej.2005.08.012), [Liu et al. (2004)](https://doi.org/10.1016/S0255-2701(02)00182-4)

| Temperature (°C) | Density (kg/m³) | Viscosity (mPa s) | Low D (10⁻⁹ m²/s) | Central D (10⁻⁹ m²/s) | High D (10⁻⁹ m²/s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 40 | 800.07 | 2.346 | 0.484 | 0.621 | 0.760 |
| 50 | 792.23 | 1.844 | 0.647 | 0.815 | 0.998 |
| 60 | 784.39 | 1.449 | 0.851 | 1.069 | 1.310 |
| 70 | 776.32 | 1.168 | 1.101 | 1.367 | 1.674 |
| 80 | 768.24 | 0.941 | 1.403 | 1.746 | 2.138 |

## Film correlation and applicability

Both correlations here use the bead diameter and the superficial velocity `u_s = Q/A` based on the empty tube cross section. The kinematic viscosity `μ/ρ` and the liquid diffusivity `D` both have units m²/s. Their ratio `Sc = μ/(ρD)` is dimensionless. Multiplying by `Re` gives the particle Péclet number `u_s d_p/D`.

Wilson and Geankoplis measured liquid mass transfer from beds of spherical benzoic-acid particles to water and aqueous propylene glycol. Their Equation 8 gives `Sh = (1.09/ε_b)(Re Sc)^(1/3)`. The paper recommends it over particle Reynolds numbers 0.0016–55 and void fractions 0.35–0.75. The authors' experiments spanned Schmidt numbers 950–70,600. The pooled comparison also includes prior data down to 165. The present bed gives `Re = 0.0117–0.1204` and `Sc ≈ 573–6060` across the five temperatures, four flows and three diffusivity scenarios. Some hot cases therefore use the part of the Schmidt range supported by the pooled literature. [Wilson and Geankoplis (1966), Equation 8, Figure 6 and nomenclature](https://doi.org/10.1021/i160017a002)

The printed reactor Equation B.2 uses `Sh = 1.17 Re^0.585 Sc^(1/3)`. With corrected `Sc`, it gives `Sh = 1.20–3.71` here, while Wilson–Geankoplis gives `Sh = 4.30–11.46`. The analytic diffusion field around one sphere in stagnant, unbounded liquid gives `Sh = 2` with diameter as the reference length. Its use is a separate low-transfer check. The source equation falls below that isolated-sphere reference in some cases, adding reason to retain the independent packed-bed correlation. [Lorenzo et al. (2025), Appendix B](https://doi.org/10.1021/acs.iecr.5c02780)
