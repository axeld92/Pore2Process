# Section-specific inert tracer rules on the saved N160 packing

Keep the same stock TRT field, geometry, balanced carrier flux, 7.032455155 g
carbon basis, matched Péclet numbers and both archived tracer time steps.
Reuse the N160 Julia particle tracker with 16,384 new inlet markers per
condition. The probe seed is 3,755,515, four archived ensemble sizes above
the original seed, so no marker seed overlaps the holdout. The four
probe conditions share this seed for paired numerical comparisons.
Sample the same flow-weighted inlet and
track continuously across seams, including recrossings. Keep the first-passage
duration between each pair of successive seams for operator construction.

Use only original sections 1–2 to initialize class mass and the first two
elapsed-time moments at the second seam, and to estimate the comparison rule
that repeats the first exchange and transit-time distribution. Retain the
original ten local transit-time ranks and their calibrated mass fractions.
Estimate each later exchange matrix and the first two moments of its receiving
section's transit time from the fresh probes. Each operator sees only the
adjacent pair of local durations. Rank independently within each section.
Use two equal, disjoint halves of the new ensemble as a sensitivity check.
Both halves keep the same archived initialization. Their spread tests
probe sampling sensitivity and does not estimate experimental uncertainty.

The local operators use the known geometry and flow within the held-out
sections. Qualification compares their composition with the independent
archived arrivals across those sections. It tests the reduced model on this
packing, with no claim of transfer to an unprobed geometry. No global age or
outlet variance from the probes enters the rules. The closure assumes earlier
elapsed time and the next transition and delay are independent within each
donor class. Propagate tracer mass and the first two arrival-time moments in
Julia. Check conservation at every seam and freeze predictions before scoring
the original sections 3–8.

The gate is a maximum absolute relative arrival-time variance error of 10%
over sections 3–8 in each size and time-step case. Report mean arrival against
the previous 5% limit and the spread between the two probe halves. Wall uptake,
adsorption and Figure 7 are outside this check.

[Spatial Markov transport](https://doi.org/10.1103/PhysRevE.78.026308)
