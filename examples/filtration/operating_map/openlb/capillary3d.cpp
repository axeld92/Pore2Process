// SPDX-License-Identifier: GPL-2.0-or-later
#include "case.h"
int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    Config c;
    for (int i=1;i<argc;++i) {
      std::string key=argv[i]; if(key=="--circle"){c.circle=true;continue;}
      if(key=="--snapshots"){c.snapshots=true;continue;}
      if(i+1==argc) throw std::invalid_argument("Missing option value.");
      std::string value=argv[++i]; std::size_t used=0;
      if(key=="--warm") {c.warm=value;continue;}
      if(key=="--resolution")c.resolution=std::stoi(value,&used);
      else if(key=="--axial")c.axial=std::stoi(value,&used);
      else if(key=="--max-steps")c.maxSteps=std::stoi(value,&used);
      else if(key=="--min-steps")c.minSteps=std::stoi(value,&used);
      else if(key=="--gradient")c.gradient=std::stod(value,&used);
      else if(key=="--tau")c.tau=std::stod(value,&used);
      else if(key=="--odd-tau")c.oddTau=std::stod(value,&used);
      else if(key=="--mu")c.mu=std::stod(value,&used);
      else if(key=="--rho")c.rho=std::stod(value,&used);
      else if(key=="--warm-scale")c.warmScale=std::stod(value,&used);
      else throw std::invalid_argument("Unknown option: "+key);
      if(used!=value.size())throw std::invalid_argument("Invalid numeric option.");
    }
    if(c.resolution<16||c.resolution>96||c.axial<4||c.axial>32||c.maxSteps<1||c.minSteps<0||c.minSteps>c.maxSteps||
       !std::isfinite(c.gradient)||c.gradient==0||!std::isfinite(c.tau)||c.tau<=.5||c.tau>1.5||
       !std::isfinite(c.oddTau)||c.oddTau<=.5||c.oddTau>10000.||
       !std::isfinite(c.mu)||c.mu<=0||!std::isfinite(c.rho)||c.rho<=0||!std::isfinite(c.warmScale))
      throw std::invalid_argument("Require resolution 16..96, axial cells 4..32, nonzero finite gradient, positive fluid properties, tau (0.5,1.5], finite warm scale and positive max steps.");
    auto mesh=createMesh(c); Model::ParametersD parameters; Model model(parameters,mesh);
    prepareGeometry(model,c); prepareLattice(model,c); setInitialValues(model,c); simulate(model,c); return 0;
  } catch(const std::exception& e) {std::cerr<<"Star channel failed: "<<e.what()<<std::endl;return 1;}
}
