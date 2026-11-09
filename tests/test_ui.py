import json
import asyncio
from zipfile import ZipFile

import panel as pn
import pytest
from bokeh.document import Document
from panel.models.vtk import VTKSynchronizedPlot

from pore2process_demo.ui import SimulationApp


def test_completed_stale_failed_zero_and_session_snapshots():
    pn.extension("vtk")
    a, b = SimulationApp(), SimulationApp()
    try:
        doc = Document()
        root = a.view.get_root(doc)
        doc.add_root(root)
        a._request_run()
        assert a.result is not None and b.result is None
        first = a.result
        assert not a.run_button.disabled and not a.download.disabled
        assert a.plotter is not None  # Real PyVista/Panel pane constructed.
        first_plotter = a.plotter
        first_vtk = list(root.select({"type":VTKSynchronizedPlot}))[0]
        assert first_vtk.scene and first_vtk.height == 420
        a.controls["radius_um"].value = 100.
        assert "previous completed run" in a.status.object
        assert a.result is first
        saved = json.loads(ZipFile(a._download()).read("summary.json"))
        assert saved["inputs"]["radius_m"] == first.inputs.radius_m
        a._request_run()
        second = a.result
        assert list(root.select({"type":VTKSynchronizedPlot}))[0] is first_vtk
        assert a.plotter is first_plotter
        assert first_vtk.scene
        assert second.Q_in_m3_s == pytest.approx(first.Q_in_m3_s*16,rel=1e-9)
        a.column_controls["length_m"].value = .04
        assert a.result is second
        assert a.column_snapshot.length_m == .04
        a.controls["radius_um"].value = 500.
        asyncio.run(a._on_run_click(None))
        assert a.result is second and "Run failed" in a.status.object
        assert not a.run_button.disabled and not a.download.disabled
        a.column_controls["length_m"].value = 0.
        assert a.download.disabled and "Column input error" in a.column_status.object
        a.column_controls["length_m"].value = .02
        a.controls["radius_um"].value = 50.
        a.controls["p_in_Pa"].value = 0.
        a._request_run()
        assert a.result.permeability_m2 is None
        assert "Unavailable" in a.metrics.object and "Unavailable" in a.column_status.object
        saved = json.loads(ZipFile(a._download()).read("summary.json"))
        assert saved["permeability_m2"] is None
        assert not a.run_button.disabled
    finally:
        a.close()
        b.close()
