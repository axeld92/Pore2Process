// SPDX-License-Identifier: GPL-2.0-or-later
// External application following the stock advectionDiffusion1d lifecycle.
// Only geometry, prescribed velocity and diagnostics are application code.
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <chrono>
using namespace olb;
using namespace olb::names;
using T=double;
using D=descriptors::D3Q27<descriptors::VELOCITY>;
template<typename V,typename DD,typename M=momenta::AdvectionDiffusionBulkTuple>
using Transport=CombinedAdvectionDiffusionRLBdynamics<V,DD,AdvectionDiffusionTRTdynamics<V,DD>,M>;
using Model=Case<AdvectionDiffusion,Lattice<T,D>>;
struct Config {int n;T side,Dm,dt,end,speed,shear,amplitude;};
T dx(const Config& c){return c.side/c.n;}
T tauOdd(const Config& c){return .5+3*c.Dm*c.dt/(dx(c)*dx(c));}
constexpr T magic=.3;
Mesh<T,3> createMesh(const Config& c){
  T h=dx(c);IndicatorCuboid3D<T> box({c.side-h,c.side-h,c.side-h},{h/2,h/2,h/2});
  Mesh<T,3> mesh(box,h,1);mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,true,true});return mesh;
}
void prepareGeometry(Model& m,const Config& c){
  auto& g=m.getGeometry();IndicatorCuboid3D<T> box({c.side,c.side,c.side},{0.,0.,0.});
  g.rename(0,1,box);g.rename(1,1);g.communicate();g.checkForErrors();g.print();
  auto& b=g.getBlockGeometry(0);
  if(b.getNx()!=c.n||b.getNy()!=c.n||b.getNz()!=c.n)throw std::runtime_error("Unexpected periodic mesh.");
}
void prepareLattice(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry();
  l.setUnitConverter<AdeUnitConverter<T,D>>(dx(c),c.dt,c.side,c.speed,c.Dm,1.);
  dynamics::set<Transport>(l,g.getMaterialIndicator(1));
  l.setParameter<descriptors::OMEGA>(1/(.5+magic/(tauOdd(c)-.5)));
  l.setParameter<collision::TRT::MAGIC>(magic);
}
void setInitialValues(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& b=m.getGeometry().getBlockGeometry(0);
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k){
    auto x=b.getPhysR({i,j,k});auto cell=l.getBlock(0).get(i,j,k);
    Vector<T,3> u{c.speed*(1+c.shear*std::cos(2*M_PI*x[1]/c.side))*c.dt/dx(c),0.,0.};
    cell.setField<descriptors::VELOCITY>(u);
    T concentration=.5+c.amplitude*std::sin(2*M_PI*x[0]/c.side);
    cell.iniEquilibrium(concentration,u.data());
  }
  l.initialize();l.setStatisticsOff();
}
struct Diagnostics {T mass=0,min=1e100,max=-1e100;int minI=0,minJ=0,minK=0,maxI=0,maxJ=0,maxK=0;};
Diagnostics diagnose(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});Diagnostics d;
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k){
    auto cell=l.getBlock(0).get(i,j,k);T value=lbm<D>::computeRho(cell);d.mass+=value;
    if(value<d.min){d.min=value;d.minI=i;d.minJ=j;d.minK=k;}
    if(value>d.max){d.max=value;d.maxI=i;d.maxJ=j;d.maxK=k;}
  }
  d.mass*=std::pow(dx(c),3);return d;
}
// The velocity and initial profile are z-invariant. Store one complete x-y
// cross-section, and check that invariance independently in the final field.
void writeField(Model& m,const Config& c,const std::string& name){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& b=m.getGeometry().getBlockGeometry(0);
  std::ofstream file(name);file<<std::setprecision(17)<<"i,j,x_m,y_m,concentration_kg_m3,ux_m_s,z_variation_kg_m3\n";
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j){
    auto cell=l.getBlock(0).get(i,j,0);T value=lbm<D>::computeRho(cell),variation=0;
    for(int k=1;k<c.n;++k){auto other=l.getBlock(0).get(i,j,k);variation=std::max(variation,std::abs(lbm<D>::computeRho(other)-value));}
    auto x=b.getPhysR({i,j,0});auto u=cell.getField<descriptors::VELOCITY>();
    file<<i<<','<<j<<','<<x[0]<<','<<x[1]<<','<<value<<','<<u[0]*dx(c)/c.dt<<','<<variation<<'\n';
  }
}
void simulate(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});std::ofstream history("inventory.csv");
  history<<std::setprecision(17)<<"step,time_s,inventory_kg,min_concentration_kg_m3,max_concentration_kg_m3,min_i,min_j,min_k,max_i,max_j,max_k\n";
  auto append=[&](int step,const Diagnostics& d){history<<step<<','<<step*c.dt<<','<<d.mass<<','<<d.min<<','<<d.max<<','<<d.minI<<','<<d.minJ<<','<<d.minK<<','<<d.maxI<<','<<d.maxJ<<','<<d.maxK<<'\n';};
  auto first=diagnose(m,c);append(0,first);writeField(m,c,"initial.csv");
  T greatest=std::max(.5-c.amplitude-first.min,first.max-(.5+c.amplitude));int peakStep=0;
  writeField(m,c,"peak.csv");const int steps=std::ceil(c.end/c.dt);auto started=std::chrono::steady_clock::now();
  for(int step=1;step<=steps;++step){
    l.collideAndStream();auto d=diagnose(m,c);append(step,d);
    T excursion=std::max(.5-c.amplitude-d.min,d.max-(.5+c.amplitude));
    if(excursion>greatest+1e-12){greatest=excursion;peakStep=step;writeField(m,c,"peak.csv");}
  }
  writeField(m,c,"final.csv");
  std::ofstream raw("raw.json");raw<<std::setprecision(17)
    <<"{\"solver\":\"OpenLB 1.9r0\",\"model\":\"stock CombinedAdvectionDiffusionRLBdynamics with AdvectionDiffusionTRTdynamics\",\"descriptor\":\"D3Q27 VELOCITY\",\"periodic\":[true,true,true],\"n\":"<<c.n
    <<",\"side_m\":"<<c.side<<",\"dx_m\":"<<dx(c)<<",\"dt_s\":"<<c.dt<<",\"steps\":"<<steps<<",\"duration_s\":"<<steps*c.dt<<",\"speed_m_s\":"<<c.speed<<",\"shear_fraction\":"<<c.shear<<",\"diffusivity_m2_s\":"<<c.Dm
    <<",\"initial_amplitude_kg_m3\":"<<c.amplitude<<",\"initial_mean_kg_m3\":0.5,\"tau_odd\":"<<tauOdd(c)<<",\"tau_even\":"<<.5+magic/(tauOdd(c)-.5)<<",\"magic\":"<<magic<<",\"peak_step\":"<<peakStep<<",\"peak_time_s\":"<<peakStep*c.dt<<",\"wall_s\":"<<std::chrono::duration<T>(std::chrono::steady_clock::now()-started).count()<<"}\n";
}
int main(int argc,char** argv){
  initialize(&argc,&argv);try{
    if(argc!=9)throw std::invalid_argument("Use N side_m D_m2_s dt_s end_s mean_velocity_m_s shear_fraction initial_amplitude_kg_m3");
    Config c{std::stoi(argv[1]),std::stod(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]),std::stod(argv[8])};
    if(c.n<16||c.n>96||!std::isfinite(c.side)||c.side<=0||!std::isfinite(c.Dm)||c.Dm<=0||!std::isfinite(c.dt)||c.dt<=0||!std::isfinite(c.end)||c.end<=0||c.end/c.dt>200000||!std::isfinite(c.speed)||c.speed<=0||!std::isfinite(c.shear)||c.shear<0||c.shear>=1||!std::isfinite(c.amplitude)||c.amplitude<0||c.amplitude>.25||c.speed*(1+c.shear)*c.dt/dx(c)>.05)
      throw std::invalid_argument("Require N16..96, positive finite SI scales, shear in [0,1), amplitude in [0,.25], max lattice velocity <=.05, <=200000 steps.");
    auto mesh=createMesh(c);Model::ParametersD p;Model model(p,mesh);
    prepareGeometry(model,c);prepareLattice(model,c);setInitialValues(model,c);simulate(model,c);return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
