# Frozen one-cell path model on a complete two-cell pulse

The [report](../../../examples/adsorption/residence_time_holdout/REPORT.md)
records cohort `20261112T044839Z-2c87f73c`. The one-cell fit was
frozen before the full two-cell comparison. The pore reference used the
same saved N96 stock TRT field, N288 tracer cells per period and 80 µs
steps. It reproduced the old early prefix and recovered 99.99998% of
the dose over the 1.7368 s window.

The path model passes the 5% arrival, main rise, mean, cumulative-profile
and fixed-time remaining-mass checks. It predicts 99% recovery 74.33 ms
late, an error of 10.19% against 10%. Its variance is 19.78% above the
reference, against 10%. Qualification fails and adsorption coupling
stays on hold. The read-only replay checked 28 sources and 42 artifacts.
