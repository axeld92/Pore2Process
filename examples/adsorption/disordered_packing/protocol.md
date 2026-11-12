# Ordered and disordered eight-cell tracer at matched carbon and Péclet numbers

Use the same box as eight axial copies of the earlier FCC cell.
The transverse extent remains one FCC cell in each periodic
direction. Keep 32 impermeable spheres of diameter 0.8 mm and
analytical void fraction 0.6 in both cases. The disordered case is
one continuous eight-section geometry, initialized from FCC and
randomized by fixed-seed periodic, nonoverlapping sphere
displacements. The sphere diameter, count, dry density and hence
the continuous carbon mass are fixed. The earlier 7.032455155 g
column inventory is the common macroscopic basis. Use the same
geometry for both particle-matched Péclet conditions, as in the
previous tracer holdout. This is one synthetic disordered packing.

Use stock OpenLB 1.9r0 D3Q27 TRT, the stock Kupershtokh force,
and the stock BounceBack and sphere indicators. The numerical
application, build and results stay outside the OpenLB tree.
For each hydraulic solve adjust the constant driving gradient to
the archived axial interstitial speed and zero mean transverse
velocity. Require a velocity-vector error below 0.2%, native
maximum Mach below 0.04 and density departure below 0.002.
Treat any higher-Mach run as an initializer and recompute its
field at a smaller time step. Check the OpenLB material mask
against independent sphere voxelization and require voxel void
fraction within 0.001 of 0.6. Record both continuous and voxel
carbon inventories. Use the existing conservative Julia face
projection with all native directional throughputs fixed.

Compare the disordered field at 48 and 64 cells per reference
length. At N64 also run a stock ordered control on one FCC cell
and repeat its field eight times. Retain the archived N96 ordered
results as a resolution control. Choose the first common positive
flow cut from the 1/16-period axial grid for the two disordered
resolutions, using only carrier fields. Begin the ordered case at
its established cell boundary. Repeat the same flux-weighted,
perfect-sink first-hit tracker and report survival, capture among
incoming survivors, initially fast survival, and fast persistence
at every reporting seam. Use 8,192 markers per new case and the
previous independent tracer inlet seed. The fast cohort is the
fastest 90% by first-section inert transit time. Use 80 and 40 µs
at N64 and 80 µs at N48. Require inert recovery of at least 99.9%.

Match each particle case's archived Péclet number
\(Pe=\bar u_\text{interstitial} L_\text{cell}/D\) using the
actual handoff speed, common reference length and a small
documented diffusivity adjustment. The measured path exchange is
allowed to respond to the geometry. Freeze the run commands and
cut before evaluating wall hits. Use the previous one-percentage-
point time-step screen on cumulative survival, conditional
capture and initially fast survival over cells 3–8. Screen the
N48 to N64 change in the same metrics at two percentage points.
If those bounds are exceeded, examine the magnitude of the
geometry effect and identify unresolved resolution sensitivity.
This milestone stops at the tracer test and uses no Figure 7
input or adsorption solve.
[OpenLB 1.9 release](https://www.openlb.net/download/) ·
[Periodicity study](https://link.springer.com/article/10.1007/s11242-021-01587-1)
