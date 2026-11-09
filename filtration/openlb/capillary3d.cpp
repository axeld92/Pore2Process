// SPDX-License-Identifier: GPL-2.0-or-later
// Application driver reusing the completed stock-component hydraulic case.
// No original OpenLB source, numerical class or geometry indicator is changed.
#include "case.h"

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    Config c;
    c.diameter=8e-6; c.length=4e-5; c.cellWidth=2e-5; c.dp=10.;
    for (int i=1;i<argc;++i) {
      std::string key=argv[i];
      if (key=="--geometry-only") { c.geometryOnly=true; continue; }
      if (i+1==argc) throw std::invalid_argument("Missing option value.");
      std::string value=argv[++i];
      std::size_t used=0;
      if (key=="--resolution") c.resolution=std::stoi(value,&used);
      else if (key=="--max-steps") c.maxSteps=std::stoi(value,&used);
      else if (key=="--radius-m") c.diameter=2*std::stod(value,&used);
      else if (key=="--length-m") c.length=std::stod(value,&used);
      else if (key=="--viscosity-Pa-s") c.mu=std::stod(value,&used);
      else if (key=="--density-kg-m3") c.rho=std::stod(value,&used);
      else if (key=="--pressure-drop-Pa") c.dp=std::stod(value,&used);
      else if (key=="--tau") c.tau=std::stod(value,&used);
      else if (key=="--steady-tolerance") c.tolerance=std::stod(value,&used);
      else throw std::invalid_argument("Unknown option: "+key);
      if (used!=value.size()) throw std::invalid_argument("Invalid complete numeric value for "+key);
    }
    if (c.resolution<8 || c.resolution>64 || c.resolution%2 || c.maxSteps<1 ||
        !std::isfinite(c.diameter) || c.diameter<8e-6 || c.diameter>1e-5 ||
        !std::isfinite(c.length) || c.length<4e-5 || c.length>8e-5 ||
        !std::isfinite(c.mu) || c.mu<=0 || !std::isfinite(c.rho) || c.rho<=0 ||
        !std::isfinite(c.dp) || c.dp==0 || std::abs(c.dp)>20 ||
        !std::isfinite(c.tau) || c.tau<=0.5 || c.tau>1.5 ||
        !std::isfinite(c.tolerance) || c.tolerance<=0 || c.tolerance>1e-3)
      throw std::invalid_argument("Require even resolution 8..64, radius 4..5 micrometres, length 40..80 micrometres, positive finite fluid properties, finite nonzero |pressure drop|<=20 Pa, tau in (0.5,1.5], tolerance in (0,1e-3], and positive max steps.");
    const T dx=c.diameter/c.resolution;
    const int nx=int(std::llround(c.length/dx));
    if (std::abs(nx*dx-c.length)>1e-9*dx || nx%4)
      throw std::invalid_argument("Length/dx must be an integer multiple of four for pressure taps.");
    auto mesh=createMesh(c);
    Model::ParametersD parameters;
    Model model(parameters,mesh);
    prepareGeometry(model,c);
    if (c.geometryOnly) {
      auto& stats=model.getGeometry().getStatistics();
      std::ofstream out("geometry.json");
      out << "{\"inlet_nodes\":" << stats.getNvoxel(3)
          << ",\"outlet_nodes\":" << stats.getNvoxel(4) << "}\n";
      return 0;
    }
    prepareLattice(model,c);
    setInitialValues(model,c);
    simulate(model,c);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Filtration support reference failed: " << e.what() << std::endl;
    return 1;
  }
}
