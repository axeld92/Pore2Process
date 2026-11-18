# Reconstruction of the published packed-bed calculation

The fitted 2025 rate law with catalyst loading per **bulk bed volume**
gets within **1.05 percentage points RMSE** of the authors' Figure 6c
curves. Applying the same volumetric rate formulation to the full 6.0 g
per **void-liquid volume** overpredicts held-out experimental conversion
by **10.34 percentage points on average**. The publication leaves the
kinetic catalyst concentration in the bed unspecified.

## Source and audit

The [source record](../../../../literature/sources/lorenzo2025/author_model.json)
has the 2025 reversible Bi–Bi Ping–Pong law inputs from Supporting Information
Eq. S1 and Table S3, the internal effectiveness factors from Table B1, and
the conditions for the 20 propanoate runs at a 2:1 alcohol-to-acid ratio.
Reactor comparisons use the digitized Figure 6c calculated lines and the
experimentally measured Figure 3c and Supporting Information Figure S1
points. Both are opened for comparison after all forward calculations.
[Lorenzo et al. (2025)](https://doi.org/10.1021/acs.iecr.5c02780)

Eq. S1 defines the positive disappearance rate of acid. With this definition,
the printed plus sign in the acid concentration equation (Eq. 2) creates
acid. Write the acid flux as `J = (u / eps_b) C - D_ax dC/dz`. Eq. 2 then
requires `dJ/dz = eta_I r_obs > 0` through the forward-reaction interval.
With Eq. 4 and `C <= C_feed`, the inlet concentration gradient is
nonpositive, so `J_in >= (u / eps_b) C_feed`. The zero outlet gradient and
depleted exit concentration require `J_out <= (u / eps_b) C_feed`. The
boundary inequalities contradict the required flux increase.
The replay explicitly uses acid consumption in Eq. 2. Eq. 4 is taken
literally in the main cases, with fixed inlet concentration and zero
outlet gradient. An inlet flux boundary checks the effect of replacing
the printed inlet by the usual Danckwerts form.

The reconstruction uses mol/L for species, g/L for immobilised catalyst
and minutes for the corresponding numeric rate constants. The publication
prints the Table S3 constants without units. The batch loading basis in
Table S2 is g/L and the minute interpretation gives a rate of 0.04056
mol/L/min at 40 °C, 25.3 g/L catalyst and fresh feed. The model assumes a
dry, product-free feed with constant liquid volume and concentrations
that change along one stoichiometric extent. The main calculation uses
reference values for the unadorned rate constants in the linear ester
term in Eq. S1. Using their temperature-dependent values alters the
result by less than 0.00005 percentage points.

Replaying the 68 post-start [batch markers](../../../../literature/sources/lorenzo2025/batch_2to1_markers.md)
gives **4.39 percentage points RMSE** at 25.3 g/L with time in minutes.
Interpreting the coefficients on an hourly basis gives 62.43 points.
The nominal 10% catalyst mass relative to the initial acid at 3.3 mol/L
implies 24.45 g/L, which gives 4.68 points. The rate law used these
batch experiments for calibration, so this checks the numerical basis
of the printed parameters. The packed-bed loading remains undetermined
by this batch comparison.

Eq. 2 divides the superficial velocity by the bed void fraction of
0.64. The geometric and void volumes are 28.1487 mL and 18.0151 mL,
consistent with the printed residence time of 18 min at 1 mL/min.
Combining Eq. 3 with Eq. 2 gives
`Pe = L (0.20 + 0.011 Re^0.48) / (eps_b^2 d_p)`. Fluid properties at
reaction conditions are unavailable. The low-Reynolds limit gives
`Pe = 115.28`, and using `Re = 1` changes conversion by at most 0.035
percentage points. The replay directly uses the printed internal
effectiveness factors of 0.852–0.884. Eq. B.2 as printed defines
`Sc = mu D / rho`, which is dimensionally inconsistent. Table B1
reports external effectiveness of 0.999, so external resistance is
omitted from the replay as in Eq. 2.

## Which loading basis explains the published curve?

For the same printed 6.0 g of immobilised catalyst, dividing by the bulk
bed volume gives 213.154 g/L and dividing by the void-liquid volume
gives 333.053 g/L. Here the rate acts on a **liquid concentration**.
Using 213.154 g/L with an interstitial residence time therefore applies
only 0.64 times the total mass exposure. In this model its exposure
equals 3.84 g on the liquid-volume basis before the effectiveness
factor. For a source formulated per bulk volume, dividing that source
by the void fraction would give the full 6.0 g exposure.

| Forward calculation | RMSE vs Figure 6c curves (pp) | RMSE vs held-out measurements (pp) |
| --- | ---: | ---: |
| Published curves read from Figure 6c | — | 4.07 |
| Bulk-volume loading with Table B1 effectiveness | 1.05 | 3.81 |
| Void-liquid loading with Table B1 effectiveness | 11.77 | 10.46 |
| Void-liquid loading with unit effectiveness | 14.83 | 13.31 |

Using the full mass exposure changes conversion by up to **13.61
percentage points** relative to the bulk-volume case. Switching to the
inlet flux boundary changes it by at most 0.27 percentage points. For
example, at 60 °C and 3 mL/min, the measured conversion is 49.7%, the
published curve gives 45.7%, the bulk-volume reconstruction gives
46.7%, and the full mass with Table B1 effectiveness gives 60.2%.

The close reconstruction supports a reduced-exposure interpretation.
The remaining up to 1.66 percentage-point deviation from the published
calculated lines exceeds the estimated plot-reading half-range of
0.5 percentage points. The authors' implementation, its exact rate units,
its bed catalyst loading basis, and its reaction-liquid transport
inputs remain uncertain. The 0.64 normalization remains an inference
about the calculation, and the physically active catalyst mass remains
unresolved.
[Comparison table](reference_comparison.csv)
· [Numerical summary](reference_summary.json)

## Reproduce

From the repository root in the configured Julia and Python environment:

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_transfer/published/test_model.jl
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transfer.published.run
```

The solver checks against a first-order batch solution and two
first-order dispersion solutions with the different inlet conditions.
The six replay scenarios use 512 and 1,024 axial cells per condition.
The maximum grid change in conversion is 0.00150 percentage points.
The worst integrated acid flux-balance error is 4.68 × 10⁻¹⁰ relative
to the inlet acid flux. Run outputs include the frozen forward grid and
the individual bed and batch comparison rows.
