# Independent liquid-film conductance from saved TRT flow

Use the saved N160 field, the eight-section periodic-side geometry and the
liquid Péclet numbers already matched to the two 25 cm operating conditions.
No breakthrough data enters the reduction. The virtual solute is removed
on its first carbon contact, with zero wall concentration and no particle
resistance. Use the flux-weighted pulse at each inlet. For a steady unit
inlet concentration, its wall capture gives uptake rate
$Q(1-S)$ and its mean first-contact-or-exit time $T_a$ gives bulk liquid
concentration $QT_a/V_l$. Their ratio gives liquid-film conductance
$K_l=V_l(1-S)/T_a$ for that domain. Use the liquid voxel volume and the
analytical external area of the same spheres to estimate $k_l=K_l/A_s$.
Use the first and second halves separately as an axial sensitivity, with
actual sphere surface areas in each half. For the second half use surviving
markers and first-passage intervals, so transport across its upstream cut
is assigned chronologically. The full-domain estimate needs no interval
assignment at internal cuts.

Use the archived 24,576 trajectories per flow at 80 µs and repeat a
matched 8,192-marker subset at 40 µs on the identical field. If the paired
point estimate changes by more than 3%, or the threshold lies within the
sampling and time-step range, continue at 20 µs. Report 95% intervals
from resampling whole trajectories and report time-step sensitivity
separately. Compare unreacted mean transit time with pore volume divided
by archived throughput. The prior N128/N160 wall-contact screen is
supporting evidence for the field, with no new spatial convergence
claim for the film coefficient.

Use the frozen 25 cm carbon mass, external sphere area based on that mass
and the same slow-access resistance for the column comparison. The
primary estimate applies the pore-scale $k_l$ to that carbon area.
Compare with a liquid-volume transfer, which applies the pore-scale
conductance per liquid volume to the same column. This alternative gauges
the consequence of using a synthetic packing with a different void
fraction and therefore a different ratio of external area to liquid
volume. Both transfers keep carbon mass unchanged. Compare with the
previous film conductance and the fixed-conductance 5% threshold, then
check the clean-film plateau with the former wall-contact allocation.
A result above the conductance threshold is necessary. The clean-film
plateau checks whether the same allocation can actually use it. No
adsorbing-bed simulations or Figure 7 fits are part of this study.
