// SPDX-License-Identifier: GPL-2.0-or-later
// Application lifecycle follows the stock poiseuille3d example.
#include "case.h"

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    Config c;
    for (int i=1;i<argc;++i) {
      std::string key=argv[i];
      if (key=="--geometry-only") { c.geometryOnly=true; continue; }
      if (i+1==argc) throw std::invalid_argument("Missing option value.");
      std::string value=argv[++i]; std::size_t used=0;
      if (key=="--resolution") c.resolution=std::stoi(value,&used);
      else if (key=="--max-steps") c.maxSteps=std::stoi(value,&used);
      else if (key=="--radius-m") c.diameter=2*std::stod(value,&used);
      else if (key=="--bulk-radius-m") c.bulkRadius=std::stod(value,&used);
      else if (key=="--lead-length-m") c.lead=std::stod(value,&used);
      else if (key=="--viscosity-Pa-s") c.mu=std::stod(value,&used);
      else if (key=="--density-kg-m3") c.rho=std::stod(value,&used);
      else if (key=="--pressure-drop-Pa") c.dp=std::stod(value,&used);
      else if (key=="--tau") c.tau=std::stod(value,&used);
      else if (key=="--steady-tolerance") c.tolerance=std::stod(value,&used);
      else throw std::invalid_argument("Unknown option: "+key);
      if (used!=value.size()) throw std::invalid_argument("Invalid complete numeric value for "+key);
    }
    if (c.resolution<8 || c.resolution>64 || c.resolution%8 || c.maxSteps<1 ||
        c.diameter!=8e-6 || !std::isfinite(c.bulkRadius) || c.bulkRadius<4e-6 || c.bulkRadius>6e-6 ||
        !std::isfinite(c.lead) || c.lead<16e-6 || c.lead>32e-6 ||
        !std::isfinite(c.mu) || c.mu<=0 || !std::isfinite(c.rho) || c.rho<=0 ||
        !std::isfinite(c.dp) || c.dp==0 || std::abs(c.dp)>20 ||
        !std::isfinite(c.tau) || c.tau<=0.5 || c.tau>1.5 ||
        !std::isfinite(c.tolerance) || c.tolerance<=0 || c.tolerance>1e-3)
      throw std::invalid_argument("Require resolution 8..64 divisible by eight, fixed 4-um neck radius, bulk radius 4..6 um, leads 16..32 um, finite positive fluid properties, nonzero |pressure|<=20 Pa, tau in (0.5,1.5], positive tolerance<=1e-3 and max steps.");
    c.length=c.bodyLength+2*c.lead;
    const T dx=c.diameter/c.resolution;
    for (T x : {c.length,c.lead,c.lead+c.bodyLength,c.span})
      if (std::abs(x/dx-std::round(x/dx))>1e-8)
        throw std::invalid_argument("Leads, taps and mesh span must lie on lattice planes.");
    auto mesh=createMesh(c);
    Model::ParametersD parameters;
    Model model(parameters,mesh);
    prepareGeometry(model,c);
    if (c.geometryOnly) {
      auto& geo=model.getGeometry();
      const int nx=int(std::llround(c.length/dx)),ny=int(std::llround(c.span/dx));
      std::ofstream voxels("voxels.csv"); voxels << "i,j,k,material\n";
      for (int i=0;i<=nx;++i) for (int j=0;j<ny;++j) for (int k=0;k<ny;++k) {
        int m=geo.get({0,i+1,j+1,k+1});
        if (m==1 || m==3 || m==4) voxels << i << ',' << j << ',' << k << ',' << m << '\n';
      }
      return 0;
    }
    prepareLattice(model,c);
    setInitialValues(model,c);
    simulate(model,c);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Graded support failed: " << e.what() << std::endl;
    return 1;
  }
}
