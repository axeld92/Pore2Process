"""Per-browser-session controls and result snapshots; synchronous Julia calls."""
from dataclasses import asdict
import html
import logging

import panel as pn

from .bridge import ColumnInputs, predict_column_flow, run_case
from .export import zip_export
from .units import from_display, microlitres_per_minute, millilitres_per_minute
from .visualization import network_pane, pressure_figure

LOG = logging.getLogger(__name__)


def display(value, unit, precision=6):
    return "Unavailable (zero network pressure drop)" if value is None else f"{value:.{precision}g} {unit}"


class SimulationApp:
    def __init__(self):
        self.result = None
        self.column_snapshot = ColumnInputs()
        self.plotter = None
        self.network_container = None
        self.controls = {
            "nx": pn.widgets.IntInput(label="Nodes along x",value=6,start=2,end=20),
            "ny": pn.widgets.IntInput(label="Nodes along y",value=4,start=1,end=8),
            "nz": pn.widgets.IntInput(label="Nodes along z",value=4,start=1,end=8),
            "spacing_mm": pn.widgets.FloatInput(label="Spacing (mm)",value=1.,step=.1),
            "radius_um": pn.widgets.FloatInput(label="Radius (µm)",value=50.,step=5.),
            "viscosity_mPa_s": pn.widgets.FloatInput(label="Viscosity (mPa s)",value=1.,step=.1),
            "p_in_Pa": pn.widgets.FloatInput(label="Inlet pressure (Pa)",value=100.,step=10.),
            "p_out_Pa": pn.widgets.FloatInput(label="Outlet pressure (Pa)",value=0.,step=10.),
        }
        self.column_controls = {
            "length_m": pn.widgets.FloatInput(label="Column length (m)",value=.02,step=.01),
            "area_m2": pn.widgets.FloatInput(label="Column area (m²)",value=1e-4,step=1e-4,format="0.000000"),
            "delta_p_Pa": pn.widgets.FloatInput(label="Column pressure drop (Pa)",value=1000.,step=100.),
        }
        self.run_button = pn.widgets.Button(label="Run simulation",color="primary",height=44,sizing_mode="stretch_width")
        self.download = pn.widgets.FileDownload(label="Download completed run (.zip)",filename="pore2process.zip",
                                                 callback=self._download_on_loop,disabled=True,sizing_mode="stretch_width")
        self.status = pn.pane.Alert("Ready. Choose inputs and run the Julia network solver.",alert_type="info",sizing_mode="stretch_width")
        self.metrics = pn.pane.Markdown("The completed run’s flow, permeability, and conservation checks will appear here.",sizing_mode="stretch_width")
        self.network_view = pn.Column(pn.pane.Markdown("Run a simulation to inspect the network."),sizing_mode="stretch_width")
        self.profile_view = pn.Column(sizing_mode="stretch_width")
        self.column_status = pn.pane.Alert("Run a network first to supply the column permeability.",alert_type="info",sizing_mode="stretch_width")
        self.run_button.on_click(self._on_run_click)
        for widget in self.controls.values():
            widget.param.watch(self._mark_stale,"value")
        for widget in self.column_controls.values():
            widget.param.watch(self._on_column_change,"value")
        self.view = self._layout()
        if pn.state.curdoc is not None:
            pn.state.on_session_destroyed(lambda context:self.close())

    def _layout(self):
        for widget in self.controls.values():
            widget.sizing_mode = "stretch_width"
            widget.width = 230
        for widget in self.column_controls.values():
            widget.width = 210
        controls = pn.Column("### Network inputs",*self.controls.values(),self.run_button,self.download,
                             width=270,styles={"background":"#f1f6f7","padding":"16px","border-radius":"12px"})
        column = pn.Column("### From pores to a Darcy column",
            "The column uses the last completed network’s permeability and viscosity. Changing these column inputs recalculates only the Julia Darcy function.",
            pn.FlexBox(*self.column_controls.values(),sizing_mode="stretch_width"),self.column_status,
            "Assumption: the same idealized repeating medium at column scale. This demonstrates upscaling; it is not a validated prediction for a real porous material.",
            sizing_mode="stretch_width")
        results = pn.Column(self.status,self.metrics,
            pn.FlexBox(pn.Column("### Pressure in the pore network",self.network_view,
                              "Drag to rotate; scroll to zoom. Coordinates and axis ticks are in mm; x is the flow direction. Tube radii are shown at their physical scale (1×).",sizing_mode="stretch_width",styles={"flex":"1 1 400px","min-width":"0"}),
                   pn.Column("### Numerical check",self.profile_view,sizing_mode="stretch_width",styles={"flex":"1 1 400px","min-width":"0"}),sizing_mode="stretch_width"),
            column,sizing_mode="stretch_width",styles={"flex":"1 1 650px","min-width":"0"})
        return pn.Column(
            pn.pane.Markdown("# Pore2Process\n**A first bridge from pore-scale flow to process-scale performance**"),
            pn.pane.Markdown("A synthetic cubic lattice of cylindrical capillaries, solved in Julia with Hagen–Poiseuille conductances. Steady, incompressible Newtonian flow; closed transverse boundaries; no gravity or junction losses. Python handles inputs, inspection, and exports. All calculations use SI units."),
            pn.FlexBox(controls,results,sizing_mode="stretch_width"),
            max_width=1600,sizing_mode="stretch_width",margin=(20,20))

    def _network_inputs(self):
        values = {name:widget.value for name,widget in self.controls.items()}
        for name, lower, upper in (("nx",2,20),("ny",1,8),("nz",1,8)):
            if not isinstance(values[name],int) or not lower <= values[name] <= upper:
                raise ValueError(f"{name} must be an integer between {lower} and {upper}.")
        return from_display(**values)

    def _column_inputs(self):
        return ColumnInputs(**{name:widget.value for name,widget in self.column_controls.items()})

    def _mark_stale(self,event=None):
        if self.result is not None:
            self.status.object = f"Inputs changed. Displayed results and download belong to previous completed run {self.result.run_id[:8]}. Press Run simulation to update."
            self.status.alert_type = "warning"

    def _request_run(self,event=None):
        if self.run_button.disabled:
            return
        try:
            inputs = self._network_inputs()
        except Exception as e:
            self._run_error(e)
            return
        self.run_button.disabled = True
        self.run_button.label = "Running…"
        self.run_button.loading = True
        self.status.object = "Solving the network in Julia. The first run also initializes and compiles the model."
        self.status.alert_type = "info"
        for widget in self.controls.values():
            widget.disabled = True
        # Yield one server tick so the browser can show the disabled button.
        # The solve still runs synchronously on this server thread, never in a worker.
        if pn.state.curdoc is not None and pn.state.curdoc.session_context is not None:
            async def complete_on_loop():
                self._execute(inputs)
            pn.state.curdoc.add_timeout_callback(complete_on_loop,50)
        else:
            self._execute(inputs)

    # Bokeh 3.10 offloads synchronous document callbacks to arbitrary workers.
    # Async callback bodies execute on its event loop. These wrappers contain
    # no await around Julia: all Julia calls remain synchronous on that thread.
    async def _on_run_click(self,event):
        self._request_run(event)

    async def _on_column_change(self,event):
        self._column_changed(event)

    async def _download_on_loop(self):
        return self._download()

    def _run_error(self,error):
        previous = "No successful result yet." if self.result is None else f"Displayed results and download still belong to completed run {self.result.run_id[:8]}."
        self.status.object = f"Run failed: {html.escape(str(error))} {previous}"
        self.status.alert_type = "danger"

    def _execute(self,inputs):
        try:
            result = run_case(inputs)
            # Build views before replacing the last successful simulation snapshot.
            figure = pn.pane.Matplotlib(pressure_figure(result),tight=True,sizing_mode="stretch_width",height=420)
            try:
                network,plotter = network_pane(result,self.network_container,self.plotter)
                self.network_container = network
            except Exception as e:
                LOG.exception("Interactive network rendering failed")
                network,plotter = pn.pane.Alert(f"3D view unavailable: {html.escape(str(e))}. Numerical results, pressure plot, and exports remain available.",alert_type="warning"),None
                self.network_container = None
            if self.plotter is not plotter:
                self.close()
            self.plotter = plotter
            self.result = result
            self.network_view.objects = [network]
            self.profile_view.objects = [figure]
            self.metrics.object = (
                f"**Completed network run `{result.run_id[:8]}`** · {len(result.pressures_Pa)} nodes · {len(result.edges)} throats\n\n"
                f"**Total inlet flow:** {display(microlitres_per_minute(result.Q_in_m3_s),'µL/min')}  \n"
                f"**Effective permeability:** {display(result.permeability_m2,'m²')}  \n"
                f"**Normalized conservation error:** {result.normalized_conservation_error:.3e}  \n"
                f"Interior imbalance: {result.max_interior_imbalance_m3_s:.3e} m³/s · Inlet/outlet mismatch: {result.flow_mismatch_m3_s:.3e} m³/s  \n"
                f"Network length: {result.length_m*1e3:g} mm · Area: {result.area_m2*1e6:g} mm²  \n"
                f"Core computation: {result.total_s*1e3:.3f} ms · Julia initialization: {result.runtime_initialization_s:.3f} s · Bridge call: {result.bridge_call_s:.3f} s")
            self.status.object = f"Completed run {result.run_id[:8]}. Results and export use this run’s saved inputs."
            self.status.alert_type = "success"
            self.download.filename = f"pore2process-{result.run_id[:8]}.zip"
            self._column_changed()
        except Exception as e:
            LOG.exception("Network run failed")
            self._run_error(e)
        finally:
            self.run_button.disabled = False
            self.run_button.loading = False
            self.run_button.label = "Run simulation"
            for widget in self.controls.values():
                widget.disabled = False

    def _column_changed(self,event=None):
        if self.result is None:
            return
        try:
            column = self._column_inputs()
            flow = predict_column_flow(self.result,column)
        except Exception as e:
            self.column_status.object = f"Column input error: {html.escape(str(e))}. Correct it to re-enable export; network result is unchanged."
            self.column_status.alert_type = "danger"
            self.download.disabled = True
            return
        self.column_snapshot = column
        self.column_status.object = (f"Column flow: {display(millilitres_per_minute(flow),'mL/min')} · "
                                     f"Permeability and viscosity from network run {self.result.run_id[:8]}.")
        self.column_status.alert_type = "info" if flow is None else "success"
        self.download.disabled = False

    def _download(self):
        if self.result is None or self.download.disabled:
            raise ValueError("Complete a valid network and column configuration before exporting.")
        return zip_export(self.result,self.column_snapshot)

    def close(self):
        if self.plotter is not None:
            self.plotter.close()
            self.plotter = None
