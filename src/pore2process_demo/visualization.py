"""Scientific plots and the PyVista geometry used by exports and Panel."""
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / ".local/matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[2] / ".local/cache"))
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import pyvista as pv


def pressure_profile(result):
    x, inverse = np.unique(result.coordinates_m[:, 0], return_inverse=True)
    mean = np.bincount(inverse, weights=result.pressures_Pa) / np.bincount(inverse)
    return x, mean


def pressure_figure(result):
    figure = Figure(figsize=(6.2, 3.6), layout="constrained")
    FigureCanvasAgg(figure)
    ax = figure.subplots()
    x, mean = pressure_profile(result)
    p = result.inputs
    analytical = p.p_in_Pa - (p.p_in_Pa-p.p_out_Pa)*x/result.length_m
    ax.plot(x*1e3, analytical, color="#ed9a35", lw=2, label="Analytical")
    ax.plot(x*1e3, mean, "o", color="#137c83", ms=6, label="Julia network mean")
    ax.set(xlabel="x (mm)", ylabel="Pressure (Pa)", title="Pressure along the sample")
    ax.grid(alpha=.2)
    ax.legend(frameon=False)
    return figure


def network_mesh(result):
    lines = np.column_stack((np.full(len(result.edges), 2, dtype=np.int64), result.edges)).ravel()
    mesh = pv.PolyData(result.coordinates_m.copy(), lines=lines)
    mesh.point_data["node_id_zero_based"] = np.arange(len(result.pressures_Pa), dtype=np.int64)
    mesh.point_data["pressure_Pa"] = result.pressures_Pa
    mesh.cell_data["flow_from_first_to_second_m3_s"] = result.flows_m3_s
    mesh.cell_data["radius_m"] = result.radii_m
    mesh.cell_data["length_m"] = result.lengths_m
    mesh.field_data["coordinate_units"] = ["m"]
    mesh.field_data["flow_orientation"] = ["Positive flow from first to second endpoint; zero-based node IDs"]
    mesh.field_data["run_id"] = [result.run_id]
    return mesh


def network_pane(result, existing=None, plotter=None):
    import panel as pn

    mesh = network_mesh(result)
    mesh.points *= 1e3
    mesh.point_data["Pressure (Pa)"] = result.pressures_Pa
    mesh.set_active_scalars("Pressure (Pa)")
    tubes = mesh.tube(radius=result.inputs.radius_m*1e3, n_sides=8)
    if plotter is None:
        plotter = pv.Plotter(off_screen=True, window_size=(720, 430))
    else:
        plotter.clear()
    plotter.set_background("#f6f9fb")
    lo, hi = float(result.pressures_Pa.min()), float(result.pressures_Pa.max())
    if lo == hi:
        lo, hi = lo-1, hi+1
    plotter.add_mesh(tubes, scalars="Pressure (Pa)", cmap="viridis", clim=(lo,hi),
                     show_scalar_bar=False, smooth_shading=False)
    plotter.view_isometric()
    plotter.reset_camera()
    axes = {"fontsize": 12, "digits": 2}
    for index, name in enumerate(("xticker", "yticker", "zticker")):
        low,high=mesh.points[:,index].min(),mesh.points[:,index].max()
        values = low+(high-low)*np.array([.25,.5,.75])
        values = np.unique(values)
        axes[name] = {"ticks": values.tolist(), "labels": [f"{'xyz'[index]} {v:.2f}" for v in values]}
    if existing is None:
        pane = pn.pane.VTK(plotter.ren_win, height=420, sizing_mode="stretch_width",
                           orientation_widget=False, axes=axes)
        container = pn.Column(pane, pane.construct_colorbars(), sizing_mode="stretch_width")
    else:
        # Keep the mounted vtk.js pane across runs. Replacing nested VTK panes
        # can leave the browser view unmounted after the previous view closes.
        container = existing
        pane = container[0]
        # Keep the render window too: replacing it can leave vtk.js resize
        # observers attached to a disposed orientation widget.
        if pane.object is plotter.ren_win:
            pane.synchronize()
        else:
            pane.object = plotter.ren_win
        pane.axes = axes
        pane.color_mappers = list(pane.get_color_mappers(infer=True))
        container[1].object = pane.construct_colorbars().object
        pane.reset_camera()
    return container, plotter
