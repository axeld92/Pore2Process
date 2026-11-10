// SPDX-License-Identifier: GPL-2.0-or-later
#include "case.h"

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    Config c;
    for (int i=1;i<argc;++i) {
      std::string key=argv[i];
      if (key=="--geometry-only") { c.geometryOnly=true; continue; }
      if (i+1==argc) throw std::invalid_argument("Missing option value.");
      std::string value=argv[++i];
      std::size_t consumed=0;
      if (key=="--resolution") c.resolution=std::stoi(value,&consumed);
      else if (key=="--max-steps") c.maxSteps=std::stoi(value,&consumed);
      else if (key=="--radius-m") c.diameter=2*std::stod(value,&consumed);
      else if (key=="--viscosity-Pa-s") c.mu=std::stod(value,&consumed);
      else if (key=="--density-kg-m3") c.rho=std::stod(value,&consumed);
      else if (key=="--pressure-drop-Pa") c.dp=std::stod(value,&consumed);
      else if (key=="--tau") c.tau=std::stod(value,&consumed);
      else if (key=="--steady-tolerance") c.tolerance=std::stod(value,&consumed);
      else throw std::invalid_argument("Unknown option: "+key);
      if (consumed!=value.size()) throw std::invalid_argument("Invalid complete numeric value for "+key);
    }
    if (c.resolution<8 || c.resolution>64 || c.resolution%2 || c.maxSteps<1 ||
        !std::isfinite(c.diameter) || c.diameter<4e-5 || c.diameter>1.2e-4 ||
        !std::isfinite(c.mu) || c.mu<=0 || !std::isfinite(c.rho) || c.rho<=0 ||
        !std::isfinite(c.dp) || c.dp==0 || std::abs(c.dp)>1 ||
        !std::isfinite(c.tau) || c.tau<=0.5 || c.tau>1.5 ||
        !std::isfinite(c.tolerance) || c.tolerance<=0 || c.tolerance>1e-3)
      throw std::invalid_argument("Require even resolution 8..64, radius 20..60 micrometres, positive finite viscosity/density, finite nonzero |pressure drop|<=1 Pa, tau in (0.5,1.5], tolerance in (0,1e-3], and positive max steps.");
    const T dx=c.diameter/c.resolution;
    const int nx=int(std::llround(c.length/dx));
    if (std::abs(nx*dx-c.length)>1e-10*dx || nx%4)
      throw std::invalid_argument("Length/dx must be an integer multiple of four for pressure taps.");
    // Follow the official poiseuille3d example lifecycle.
    auto mesh=createMesh(c);
    Model::ParametersD parameters;
    Model model(parameters,mesh);
    prepareGeometry(model,c);
    if (c.geometryOnly) {
      auto& stats=model.getGeometry().getStatistics();
      std::ofstream out("geometry.json");
      out << "{\"resolution\":" << c.resolution << ",\"expected_boundary_nodes\":" << stats.getNvoxel(3)
          << ",\"inlet_nodes\":" << stats.getNvoxel(3) << ",\"outlet_nodes\":" << stats.getNvoxel(4) << "}\n";
      return 0;
    }
    prepareLattice(model,c);
    setInitialValues(model,c);
    simulate(model,c);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Hydraulic benchmark failed: " << error.what() << std::endl;
    return 1;
  }
}
