# Van Hecke et al. (2024): independently measured acid-inhibition shape

[Open-access original paper](https://doi.org/10.1021/acs.oprd.4c00274),
*Organic Process Research & Development* 28, 3989–4002.

Section 3.3, Figures 6–7 and Tables 4, 8 investigate solvent-free
propionic acid + geraniol with Novozym 435 at 60 °C. Four forward
progress profiles use sieves to remove water; three additional
equilibrium runs omit them. Section 3.5 rewrites the empirical
inhibition denominator so that its half-rate acid concentration
is in mol/kg. For propionic acid and geraniol, Table 8 gives
exponent 3.5 and half-rate concentration 1.70 mol/kg.
[`acid_inhibition.json`](acid_inhibition.json) records the target
comparison. These are parameters estimated from reaction profiles,
not directly measured acid dissociation or intrinsic enzyme binding
constants. The study used low catalyst loadings and a different
alcohol. The source supports a strong, reversible acid-dependent
rate effect and suggests prolonged acid exposure may cause
irreversible loss; it does not establish a dynamic deactivation
constant for pentyl propanoate or its packed bed.

In the new target model, acid concentration is calculated as
moles of acid per kilogram of conserved dry-feed liquid mass.
Use the published concentration as a fixed analogue sensitivity,
and fit the same shape with independent target batch data as a
separate candidate. Rate amplitude and apparent equilibrium
always come from the target training data. Do not label the
fitted coefficient as a transferable fundamental Ki. A lumped
product response cannot distinguish water from ester along a
dry-feed stoichiometric reaction trajectory.
