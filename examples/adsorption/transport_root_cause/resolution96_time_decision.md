# Interpreting the conditional N96 time-step test

This decision was recorded after the session restart, before any N96 scalar
execution or result. It clarifies the prospective [N96 extension](resolution96_protocol.md)
without changing the original simulation or engineering tolerance.

Run the matched half-time-step replay if the N96 largest excursion passes the
0.02 kg/m³ engineering screen or falls to at most half the matched N64
excursion of 0.1285433555070885 kg/m³. The latter is a trigger for more
verification, not acceptance of an excursion above 0.02 kg/m³.

If the finer case still has a larger excursion than both triggers, retain its
actual peak field, first two local population audits and final mixing measures.
The existing N64 half-time-step trial establishes temporal sensitivity in this
configuration. No new N96 time test is needed to report failure of the tested
matched-time-step spatial series. No periodic case alone qualifies a column
dispersion, outlet curve or adsorption prediction.
