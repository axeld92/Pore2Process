# Browser verification — 7 November 2026

Checked the locally served Panel app at `http://127.0.0.1:5006/app` with the
available computer-use tool. Initial rendering was observed in Codex's in-app
browser. After that browser's reload/control calls timed out, interaction checks
continued in a new Firefox tab, preserving the existing tab. HTTP status 200 was
also checked separately; it was not counted as browser rendering evidence.

| Check | Evidence / result |
| --- | --- |
| Default Run | 96 nodes, 224 throats, 47.1239 µL/min, 2.45437e-12 m², default column 0.736311 mL/min |
| Busy state | Run and network controls disabled while solving, restored on success/failure |
| Input edits | Changing inlet pressure to 200 Pa displayed a warning identifying the previous completed run; no automatic solve |
| Second Run | 94.2478 µL/min at 200 Pa, unchanged permeability; repeated runs complete after the event-loop callback fix |
| 3D rendering | Pressure-colored physical-radius tubes, readable axis ticks, orientation widget and labelled Pressure (Pa) colorbar |
| Repeated 3D update | The mounted VTK pane is reused; network stays visible and the colorbar updates to 0–200 Pa |
| Interactive zoom | Wheel scrolling over the VTK view visibly changed the network's scale |
| Pressure profile | Numerical points coincide with the analytical line; changed-pressure profile ranges from 200 to 0 Pa, with mm/Pa axes |
| Invalid network | 500 µm radius at 1 mm spacing rejected because diameter must be smaller than spacing; last successful run/results remain identified and Run recovers |
| Zero drop | 0/0 Pa produces zero flow and zero conservation diagnostics; inferred permeability and column flow say Unavailable |
| Column-only edit | .02 to .04 m halves column flow, keeping network run ID unchanged |
| Browser download | Saved `results/pore2process-cd92d4ad.zip`; independent ZIP/JSON read verified all five files, 200 Pa inlet, .04 m column, matching flow and network run ID |
| Input snapshot | Browser archive describes completed inputs; Python tests additionally verify that stale controls cannot overwrite that snapshot |

Screenshots `browser_changed.png` and `browser_zero.png` show the changed and
zero-flow states. `browser_verified.png` shows the final default state.
They are cropped from actual browser screenshots to the app and address bar,
excluding unused desktop space. No simulated browser image was used.

The live checks exposed and resolved two issues: Bokeh 3.10's automatic worker
dispatch for synchronous callbacks (Julia calls now remain on the main thread),
and disappearing nested VTK views when replacing their mounted pane (the pane is
now reused). A Bokeh Document-backed regression test covers repeated VTK updates.
Validation messages omit Julia stack traces from the interface; detailed traces
remain in the local server log.

The original native automation's drag attempts selected surrounding text and
did not conclusively exercise drag-to-rotate. This was closed on 9 November:
an actual canvas drag in Codex's in-app browser visibly rotated the network.
`browser_capillary_rotated.jpg` records that view. The original capillary
default was rerun and showed 47.1239 µL/min and 2.45437e-12 m². Interactive
rendering, zoom, changes of scientific inputs, zero/invalid states, and browser
export were checked as described above.
No mobile viewport, multi-user load, offline JS-asset installation, other OS,
or production deployment was tested. Session-state independence is covered by
the Python UI test, rather than a browser load test.

## Explicit-channel and filtration interface — 9 November 2026

The final app was checked in Codex's in-app browser and Firefox at
`http://127.0.0.1:5006/app`. The accepted native evidence opens the filtration
tab by default; the original capillary reference remains in the other tab.
No HTTP/import check is counted as a visual pass.

| Check | Observed result |
| --- | --- |
| Default filtration run | Actual Julia translation and DWSIM child solve/save/reload; 12 LMH, inlet Re=20, mean TMP 1.221 kPa, shaft energy 0.000912569 kWh/m³ |
| Busy controls | Run and inputs disabled during execution and restored on completion; no automatic solve on edits |
| Changed support | Selecting 6 µm marks the previous result stale; a new actual run gives 1.083 kPa, 0.000810323 kWh/m³ and 11.20% saving against the verified matching straight-support case |
| Profiles | Labelled TMP/flux/COD plots update to the completed case; sustainable flux remains explicitly unavailable |
| Native 3D view | Actual low-gradient OpenLB velocity field, readable x/y/z ticks in mm, labelled m/s colour bar and physical-scale geometry; an actual canvas drag visibly rotates it |
| Invalid input | Empty flux becomes an invalid numeric input; Run reports `flux_LMH must be finite and positive`, restores controls and retains the prior completed result/download. The 8–12 LMH widget bounds also constrain ordinary out-of-range entry |
| Browser isolation | A new browser session starts with default 4 µm controls and no process result while another session retains its completed 6 µm case |
| Browser download | Firefox's native Save dialog saved `results/pore2process-filtration-20261109T164208Z-9e81f7d0.zip`: 363,257 bytes, 11 files, all manifest hashes correct |
| Download contents | Correct 6 µm / 12 LMH snapshot/run ID, JSON nulls for unavailable quantities, 101 SI profile rows, 1,860 native VTP points with correct x=dx and zero-based IDs, actual DWSIM archive and labelled plots |
| Failed-run download | The complete link remains enabled with the same blob URL and completed run after invalid input; no partial or NaN input snapshot replaces it |
| Original reference | Default 47.1239 µL/min and 200 Pa case 94.2478 µL/min with unchanged k; pressure colour bar updates to 0–200 Pa and the network remains visible |
| Two VTK views | Both stay mounted during tab switches; repeated network runs preserve the render window. The completed process remains available after switching back |
| Diagnostics | No browser JavaScript errors during the final checks; no server traceback. A Bokeh colour-bar-only plot emits its missing-renderers startup warning |

Live checks exposed a VTK unmount/optional-orientation-marker resize issue and
a stale empty download link when an embedded FileDownload was updated while
disabled. The final app keeps its VTK panes/render window, uses explicit tick
labels without the optional orientation marker, and mounts a new complete
embedded download widget only after each successful archive is ready. No vendor
JavaScript is patched. The failed zero-byte probe is retained under
`results/star_channel/diagnostics/`; it is not the accepted download.

In-app download-event capture timed out, so download acceptance uses the actual
Firefox Save dialog and independent file inspection. Screenshots
`browser_filtration_verified.jpg`, `browser_filtration_3d.jpg` and
`browser_capillary_tabs.jpg` show actual final views. The
[browser provenance record](star_channel_browser_reference.json) retains the
download SHA, final UI source hashes, screenshot hashes and separate HTTP 200
observation. UI presentation fixes followed the immutable numerical matrix;
their final source hashes are recorded separately. The final Python suite
passes 61 tests; Julia's 878 assertions and the accepted native/process checks
are recorded in [the scientific report](star_channel_showcase.md).

No mobile viewport, multi-user load, offline JavaScript asset installation,
other OS or production deployment is claimed. The biological retention remains
a scenario assumption and no sustainable-flux/fouling-time validation is claimed.
