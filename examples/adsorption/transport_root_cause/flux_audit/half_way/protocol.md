# Prospective finite stock half-way reflection candidate

Change only the scalar wall realization: retain stock D3Q27, CombinedRLB/TRT,
FirstOrder equilibrium, physical D=1e-9 m2/s, velocity, dx/dt and elapsed window.
Use stock `setBouzidiBoundary<T,D,BouzidiPostProcessor>` to register fluid-solid
link postprocessors, using stock analytical geometry indicators. Configure only
the already active BOUZIDI_DISTANCE entries to exactly .5; keep inactive entries
at -1. This disables curved-distance interpolation and retains the same voxel
midlink wall as full-way BounceBack, while removing its extra return-step lag.
Solids use stock NoDynamics. Do not implement a boundary/collision operator.

After l.initialize, restore declared fluid initial populations using stock
iniEquilibrium(C,u); initialize invokes postprocessing and can otherwise alter
the initial field. Require complete fluid-to-solid link coverage, no reflection
on fluid-fluid links, exact q=.5, untouched inactive entries and periodic seams.
The stock PostStream operator returns a copied outgoing population immediately
to fluid; its copy on a solid node is NOT additional inventory. Require physical
fluid mass conservation and exact initial concentration snapshots. Read-only
diagnostics may compare actual reflected populations and local fluid balance,
accounting explicitly for the intentional PostStream update.

Stage1: flat-wall constant32, diffusion-only32, Poiseuille32 and Poiseuille96,
four z layers, with the previously qualified4096-point Neumann continuum oracle.
Keep that benchmark's H/U/dt/11294steps and amplitude. Preserve every-step fluid
inventory/extrema and initial/peak/final fields. Structural, finite, mass,
constant, diffusion and N96 profile/bounds screens must pass before Stage2.
N32 Poiseuille profile failure is an expected unresolved coarse-grid control,
not a reason to relax screens. N96 must retain <=1% RMS/amplitude.

Stage2, only if Stage1 qualifies: periodic FCC sine32 and sine64 with the
already accepted, same-grid native velocity fields and original matched dt/
window. No pressure/velocity fit, changed diffusivity, clipping, normalization,
dispersion coefficient or process prediction. Preserve bound failures. At most
six completed native candidate cases; do not expand to another algorithm if
this finite candidate fails. Use separate project code/build/results folders,
source/artifact/executable hashes, actual commands and unchanged-vendor checks.
Independent review of stock API/lifecycle, coverage, inventory and reference
is performed by the channel-control agent. No OpenLB vendor file is modified.
