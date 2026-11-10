// SPDX-License-Identifier: GPL-2.0-or-later
// External application using OpenLB 1.9r0 advectionDiffusion1d's lifecycle.
// The filename reuses the existing external build helper; this is a 2D ADE test.
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <chrono>
using namespace olb;
using namespace olb::names;
using T=double;
using D=descriptors::D2Q5<descriptors::VELOCITY>;
using Model=Case<AdvectionDiffusion,Lattice<T,D>>;
struct Config { int n=32; T pe=4.,tauOdd=.8,amplitude=.2; };
constexpr T length=1e-3, diffusivity=9.1e-10, magic=.1875;

Mesh<T,2> createMesh(const Config& c) {
  const T h=length/c.n;
  IndicatorCuboid2D<T> domain({length-h,3*h},{0.,0.});
  Mesh<T,2> mesh(domain,h,1);mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,true});return mesh;
}
void prepareGeometry(Model& model) {
  auto& g=model.getGeometry();g.rename(0,1);g.communicate();g.checkForErrors();g.print();
}
void prepareLattice(Model& model,const Config& c) {
  auto& l=model.getLattice(AdvectionDiffusion{});auto& g=model.getGeometry();
  const T h=length/c.n,dt=(c.tauOdd-.5)*h*h/(3*diffusivity);
  l.setUnitConverter<AdeUnitConverter<T,D>>(h,dt,length,c.pe*diffusivity/length,diffusivity,1.);
  dynamics::set<AdvectionDiffusionTRTdynamics>(l,g,1);
  // Scalar diffusivity is controlled by the ODD relaxation channel. Stock TRT's
  // OMEGA parameter is the EVEN rate; do not pass the BGK converter rate here.
  const T tauEven=.5+magic/(c.tauOdd-.5);
  l.setParameter<descriptors::OMEGA>(1/tauEven);
  l.setParameter<collision::TRT::MAGIC>(magic);
}
void setInitialValues(Model& model,const Config& c) {
  auto& l=model.getLattice(AdvectionDiffusion{});auto& g=model.getGeometry();
  const auto& cv=l.getUnitConverter();
  AnalyticalConst2D<T,T> velocity(cv.getCharLatticeVelocity(),0.);
  fields::set<descriptors::VELOCITY>(l,g.getMaterialIndicator({1}),velocity);
  auto& b=g.getBlockGeometry(0);T u[2]={cv.getCharLatticeVelocity(),0.};
  if (b.getNx()!=c.n) throw std::runtime_error("Unexpected periodic grid length.");
  for(int i=0;i<b.getNx();++i)for(int j=0;j<b.getNy();++j){
    auto x=b.getPhysR({i,j});
    l.getBlock(0).get(i,j).iniEquilibrium(1+c.amplitude*std::sin(2*M_PI*x[0]/length),u);
  }
  l.initialize();
}
void getResults(Model& model,const Config& c,int steps,T wall) {
  auto& l=model.getLattice(AdvectionDiffusion{});auto& b=model.getGeometry().getBlockGeometry(0);
  const auto& cv=l.getUnitConverter();
  l.setProcessingContext(ProcessingContext::Evaluation);
  SuperLatticeDensity2D<T,D> concentration(l);
  std::ofstream out("field.csv");out<<std::setprecision(17)<<"i,j,x_m,y_m,normalized_concentration\n";
  for(int i=0;i<b.getNx();++i)for(int j=0;j<b.getNy();++j){
    const int pos[3]={0,i,j};T rho[1];if(!concentration(rho,pos))throw std::runtime_error("Field read failed.");
    auto x=b.getPhysR({i,j});out<<i<<','<<j<<','<<x[0]<<','<<x[1]<<','<<rho[0]<<'\n';
  }
  std::ofstream raw("raw.json");raw<<std::setprecision(17)
    <<"{\"solver\":\"OpenLB 1.9r0\",\"dynamics\":\"stock AdvectionDiffusionTRTdynamics\",\"descriptor\":\"D2Q5\",\"nx\":"<<b.getNx()<<",\"ny\":"<<b.getNy()
    <<",\"length_m\":"<<length<<",\"diffusivity_m2_s\":"<<diffusivity<<",\"velocity_m_s\":"<<cv.getCharPhysVelocity()
    <<",\"dx_m\":"<<cv.getPhysDeltaX()<<",\"dt_s\":"<<cv.getPhysDeltaT()<<",\"steps\":"<<steps
    <<",\"tau_odd\":"<<c.tauOdd<<",\"tau_even\":"<<.5+magic/(c.tauOdd-.5)<<",\"magic\":"<<magic
    <<",\"initial_amplitude\":"<<c.amplitude<<",\"wall_s\":"<<wall<<"}\n";
}
void simulate(Model& model,const Config& c) {
  auto& l=model.getLattice(AdvectionDiffusion{});const auto& cv=l.getUnitConverter();
  const int steps=std::lround(.05*length*length/diffusivity/cv.getPhysDeltaT());
  auto start=std::chrono::steady_clock::now();
  for(int i=0;i<steps;++i)l.collideAndStream();
  getResults(model,c,steps,std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count());
}
int main(int argc,char** argv){
  initialize(&argc,&argv);
  try{
    if(argc!=5)throw std::invalid_argument("Use: resolution Peclet tau_odd amplitude");
    Config c{std::stoi(argv[1]),std::stod(argv[2]),std::stod(argv[3]),std::stod(argv[4])};
    if(c.n<16||c.n>256||!std::isfinite(c.pe)||std::abs(c.pe)>10||!std::isfinite(c.tauOdd)||c.tauOdd<=.5||c.tauOdd>1.5||!std::isfinite(c.amplitude)||c.amplitude<0||c.amplitude>.5)
      throw std::invalid_argument("Require N=16..256, |Pe|<=10, tau_odd in (0.5,1.5], amplitude=0..0.5.");
    auto mesh=createMesh(c);Model::ParametersD p;Model model(p,mesh);
    prepareGeometry(model);prepareLattice(model,c);setInitialValues(model,c);simulate(model,c);return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
