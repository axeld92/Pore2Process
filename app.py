"""Panel executes this file separately for each browser session."""
import panel as pn
from pore2process_demo.ui import SimulationApp
from pore2process_demo.filtration_ui import FiltrationApp

pn.extension("vtk")
simulation = SimulationApp()
filtration = FiltrationApp()
pn.Tabs(('Capillary network',simulation.view),('Filtration showcase',filtration.view),
        active=1 if filtration.ready else 0,dynamic=False,
        sizing_mode='stretch_width').servable(title="Pore2Process · Pores to process")
