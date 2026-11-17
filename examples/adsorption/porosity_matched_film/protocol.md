# Clean-wall transfer at 42.7% void fraction

Keep the 25 cm column's 17.581 g of 0.80 mm carbon and its frozen
slow-access resistance. Make a separate synthetic packing of 32
nonoverlapping 0.80 mm spheres with periodic sides and analytical
void fraction 0.427. Preserve approximately the eight-section sample
length while reducing the transverse area. Seed one randomized geometry
and keep the same centers at every flow and grid. A single periodic
packing is a scenario for this void fraction.

Use the existing external stock OpenLB 1.9r0 application with D3Q27,
TRT, Kupershtokh forcing and BounceBack. Adapt the dimensions and fluid
viscosity in a separate generated source directory. Use separate native
flows at the two column interstitial velocities. The native checks
retain 0.2% velocity, 0.04 maximum Mach and 0.002 maximum relative density
departure limits. Verify the voxel inventory against the analytical
geometry. Balance the saved faces with the existing Julia projection
and require its existing flux and conservation limits.

Match axial-section Péclet numbers 13,291.35 and 35,652.76 by deriving
the archive tracer diffusivity from measured throughput and void volume
on each grid. Use the same cut across the compared cases when possible.
Check recovery of inert markers, their mean residence against liquid
volume divided by flow, and clean-wall first contact at 40 and 20 µs.
Compare transverse N45 and N60 at the fine tracer step. The time and
grid screens are a 3% relative change in film conductance and a
2-percentage-point absolute change in projected escape. Use finer runs
if the screens miss and the 5% outcome still depends on resolution.
Use 6,144 markers per flow on N45 and 8,190 per flow and step on N60.
For the path model, reserve the first third for the initial transit-time
classification and the remaining two thirds for exchange and wall contact.

Apply the first-contact uptake divided by mean liquid concentration
method from the previous film test. Scale the film by external carbon
area and by liquid volume as separate checks, retain the fixed slow
resistance, and recalculate the wall-contact distribution from the new
flow. Report projected 25 cm clean-bed escape and tracer sampling
intervals for both rates. Count a clear result below 5% at high flow
only if the relevant sampling interval is below the cutoff and the
numerical comparison supports it. Report any unresolved sensitivity.
No carbon loading or Figure 7 fit enters the test.
