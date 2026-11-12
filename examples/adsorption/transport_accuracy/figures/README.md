# Saved-data figures

`uniform_accuracy.png` distinguishes local extrema from average physical-liquid
concentration through sixteen void-volume transits. `front_diagnostics.png`
compares native concentration bounds and absolute outlet flux for the sharp
input, one-transit gradual input and stock regularized-inlet diagnosis. Upper
and lower concentration envelopes are shown; no values are clipped.
`periodic_isolation.png` compares the actual free-flow field with its exact
solution and shows the porous extrema against the initial .25–.75 bounds.
`bounded_diagnosis.png` shows the four subsequent controls and the reused
reference, including their largest excursions beyond that range. Its separate
hash and inspection record is `results/adsorption/transport_accuracy/diagnosis_figure_sources.json`.

Plots use Matplotlib and saved native histories. Exact inputs/PNG hashes and
visual inspection status are recorded in
`results/adsorption/transport_accuracy/figure_sources.json`. These are scientific
figures, not browser/interface acceptance.
