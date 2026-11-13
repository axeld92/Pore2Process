# Two-section disordered exchange holdout

Cohort `20261113T170400Z-e462ad77` used the first two N160
sections to estimate inert class exchange and wall-free transfer
conditional on the class transition. Frozen Julia predictions
were then checked against sections 3–8. The wall metrics passed
in three of four runs. The 0.45 mm case at 40 µs exceeded the
fast-cohort wall-free limit at 2.14 percentage points. All
mean arrival checks passed. Variance errors of 10.47–22.06%
exceeded the 10% screen. Read-only replay verified 27 sources
and 92 artifacts. [Report and reproduction](../../../examples/adsorption/disordered_exchange/REPORT.md)

An exact variance budget attributes the section-7 error primarily to
path exchange and local transit times varying along the bed. The net
travel-time memory contribution is small there and offsets part of
the exchange effect at the outlet. The frozen prediction is unchanged.
[Variance audit](../../../examples/adsorption/disordered_exchange/VARIANCE_AUDIT.md)
