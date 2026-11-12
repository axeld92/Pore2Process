// SPDX-License-Identifier: GPL-2.0-or-later
// Flat-wall diagnostic following the stock advectionDiffusion1d lifecycle.
// Geometry, prescribed velocity and read-only accounting are application code.
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
struct Config{int n,layers,steps;T H,Dm,dt,speed,amplitude;};
constexpr T magic=.3;
T dx(const Config& c){return c.H/c.n;}
T tauOdd(const Config& c){return .5+3*c.Dm*c.dt/(dx(c)*dx(c));}
Mesh<T,3> createMesh(const Config& c){
  T h=dx(c);IndicatorCuboid3D<T> box({c.H-h,(c.n+1)*h,(c.layers-1)*h},{h/2,-h/2,h/2});
  Mesh<T,3> mesh(box,h,1);mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,false,true});return mesh;
}
void prepareGeometry(Model& m,const Config& c){
  auto& g=m.getGeometry();T h=dx(c);
  IndicatorCuboid3D<T> box({c.H,c.H+2*h,c.layers*h},{0.,-h,0.});g.rename(0,1,box);
  IndicatorCuboid3D<T> lower({c.H,.75*h,c.layers*h},{0.,-h,0.});g.rename(1,2,lower);
  IndicatorCuboid3D<T> upper({c.H,.75*h,c.layers*h},{0.,c.H+.25*h,0.});g.rename(1,2,upper);
  g.rename(1,1);g.communicate();g.checkForErrors();g.print();auto& b=g.getBlockGeometry(0);
  if(b.getNx()!=c.n||b.getNy()!=c.n+2||b.getNz()!=c.layers)throw std::runtime_error("Unexpected flat-channel mesh.");
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n+2;++j)for(int k=0;k<c.layers;++k)
    if(b.get({i,j,k})!=((j==0||j==c.n+1)?2:1))throw std::runtime_error("Flat-channel wall mask mismatch.");
}
void prepareLattice(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry();
  l.setUnitConverter<AdeUnitConverter<T,D>>(dx(c),c.dt,c.H,std::max(c.speed,T{1e-4}),c.Dm,1.);
  dynamics::set<Transport>(l,g,1);boundary::set<boundary::BounceBack>(l,g,2);
  l.setParameter<descriptors::OMEGA>(1/(.5+magic/(tauOdd(c)-.5)));
  l.setParameter<collision::TRT::MAGIC>(magic);
}
void setInitialValues(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& b=m.getGeometry().getBlockGeometry(0);
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n+2;++j)for(int k=0;k<c.layers;++k){
    auto x=b.getPhysR({i,j,k});auto cell=l.getBlock(0).get(i,j,k);T s=x[1]/c.H;
    T ux=(j>0&&j<c.n+1)?6*c.speed*s*(1-s):0.;Vector<T,3> u{ux*c.dt/dx(c),0.,0.};
    cell.setField<descriptors::VELOCITY>(u);T value=.5+c.amplitude*std::sin(2*M_PI*x[0]/c.H);cell.iniEquilibrium(value,u.data());
  }
  l.initialize();l.setStatisticsOff();
}
struct Diagnostics{T fluid=0,wall=0,min=1e100,max=-1e100;int minI=0,minJ=0,maxI=0,maxJ=0;};
Diagnostics diagnose(Model& m,const Config& c){
  auto& block=m.getLattice(AdvectionDiffusion{}).getBlock(0);Diagnostics d;
  for(int i=0;i<c.n;++i)for(int j=1;j<=c.n;++j)for(int k=0;k<c.layers;++k){
    auto cell=block.get(i,j,k);T value=lbm<D>::computeRho(cell);d.fluid+=value;
    if(value<d.min){d.min=value;d.minI=i;d.minJ=j-1;}
    if(value>d.max){d.max=value;d.maxI=i;d.maxJ=j-1;}
  }
  // After streaming, these wall populations have arrived from physical fluid.
  // Full-way BounceBack reflects them into fluid on the following step.
  // Populations facing unused exterior nodes are not part of this closed set.
  for(int i=0;i<c.n;++i)for(int k=0;k<c.layers;++k)for(int p=0;p<D::q;++p){
    if(descriptors::c<D>(p,1)==-1){auto cell=block.get(i,0,k);d.wall+=cell[p]+descriptors::t<T,D>(p);}
    if(descriptors::c<D>(p,1)==1){auto cell=block.get(i,c.n+1,k);d.wall+=cell[p]+descriptors::t<T,D>(p);}
  }
  T volume=std::pow(dx(c),3);d.fluid*=volume;d.wall*=volume;return d;
}
std::vector<T> captureField(Model& m,const Config& c){
  auto& block=m.getLattice(AdvectionDiffusion{}).getBlock(0);std::vector<T> values(2*c.n*c.n);
  for(int i=0;i<c.n;++i)for(int j=1;j<=c.n;++j){
    auto cell=block.get(i,j,0);T value=lbm<D>::computeRho(cell),variation=0;
    for(int k=1;k<c.layers;++k){auto other=block.get(i,j,k);variation=std::max(variation,std::abs(lbm<D>::computeRho(other)-value));}
    int index=2*(i*c.n+j-1);values[index]=value;values[index+1]=variation;
  }return values;
}
void writeField(Model& m,const Config& c,const std::string& name,const std::vector<T>& values){
  auto& block=m.getLattice(AdvectionDiffusion{}).getBlock(0);auto& b=m.getGeometry().getBlockGeometry(0);
  std::ofstream out(name);out<<std::setprecision(17)<<"i,j,x_m,y_m,concentration_kg_m3,ux_m_s,z_variation_kg_m3\n";
  for(int i=0;i<c.n;++i)for(int j=1;j<=c.n;++j){
    auto cell=block.get(i,j,0);auto x=b.getPhysR({i,j,0});auto u=cell.getField<descriptors::VELOCITY>();int index=2*(i*c.n+j-1);
    out<<i<<','<<j-1<<','<<x[0]<<','<<x[1]<<','<<values[index]<<','<<u[0]*dx(c)/c.dt<<','<<values[index+1]<<'\n';
  }
}
void simulate(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});std::ofstream history("inventory.csv");history<<std::setprecision(17);
  history<<"step,time_s,fluid_inventory_kg,active_wall_population_kg,kinetic_inventory_kg,min_concentration_kg_m3,max_concentration_kg_m3,min_i,min_j,max_i,max_j\n";
  auto append=[&](int s,const Diagnostics& d){history<<s<<','<<s*c.dt<<','<<d.fluid<<','<<d.wall<<','<<d.fluid+d.wall<<','<<d.min<<','<<d.max<<','<<d.minI<<','<<d.minJ<<','<<d.maxI<<','<<d.maxJ<<'\n';};
  auto first=diagnose(m,c);append(0,first);auto peak=captureField(m,c);writeField(m,c,"initial.csv",peak);
  T greatest=std::max(.5-c.amplitude-first.min,first.max-.5-c.amplitude);int peakStep=0;auto start=std::chrono::steady_clock::now();
  for(int s=1;s<=c.steps;++s){
    l.collideAndStream();auto d=diagnose(m,c);append(s,d);T error=std::max(.5-c.amplitude-d.min,d.max-.5-c.amplitude);
    if(error>greatest+1e-12){greatest=error;peakStep=s;peak=captureField(m,c);}
  }
  writeField(m,c,"peak.csv",peak);writeField(m,c,"final.csv",captureField(m,c));
  std::ofstream raw("raw.json");raw<<std::setprecision(17)<<"{\"solver\":\"OpenLB 1.9r0\",\"model\":\"stock CombinedRLB with stock TRT, FirstOrder equilibrium, fullway BounceBack\",\"descriptor\":\"D3Q27 VELOCITY\",\"periodic\":[true,false,true],\"n\":"<<c.n<<",\"z_layers\":"<<c.layers<<",\"height_m\":"<<c.H<<",\"length_m\":"<<c.H<<",\"depth_m\":"<<c.layers*dx(c)<<",\"fluid_volume_m3\":"<<c.H*c.H*c.layers*dx(c)<<",\"wall_locations_m\":[0,"<<c.H<<"],\"dx_m\":"<<dx(c)<<",\"dt_s\":"<<c.dt<<",\"steps\":"<<c.steps<<",\"duration_s\":"<<c.steps*c.dt<<",\"mean_speed_m_s\":"<<c.speed<<",\"converter_characteristic_speed_m_s\":"<<std::max(c.speed,T{1e-4})<<",\"sampled_maximum_lattice_speed\":"<<1.5*c.speed*(1-1./(c.n*c.n))*c.dt/dx(c)<<",\"diffusivity_m2_s\":"<<c.Dm<<",\"initial_mean_kg_m3\":0.5,\"initial_amplitude_kg_m3\":"<<c.amplitude<<",\"tau_odd\":"<<tauOdd(c)<<",\"tau_even\":"<<.5+magic/(tauOdd(c)-.5)<<",\"magic\":"<<magic<<",\"peak_step\":"<<peakStep<<",\"peak_time_s\":"<<peakStep*c.dt<<",\"maximum_excursion_kg_m3\":"<<greatest<<",\"initial_active_wall_population_kg\":"<<first.wall<<",\"wall_s\":"<<std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}
int main(int argc,char** argv){
  initialize(&argc,&argv);try{
    if(argc!=9)throw std::invalid_argument("Use N z_layers H_m D_m2_s dt_s steps mean_speed_m_s initial_amplitude_kg_m3");
    Config c{std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[6]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[7]),std::stod(argv[8])};
    if(c.n<16||c.n>256||c.layers<4||c.layers>96||c.steps<1||c.steps>200000||!std::isfinite(c.H)||c.H<=0||!std::isfinite(c.Dm)||c.Dm<=0||!std::isfinite(c.dt)||c.dt<=0||!std::isfinite(c.speed)||c.speed<0||!std::isfinite(c.amplitude)||c.amplitude<0||c.amplitude>.25||1.5*c.speed*c.dt/dx(c)>.05)
      throw std::invalid_argument("Require N16..256, z_layers4..96, positive finite H/D/dt, nonnegative U, amplitude[0,.25], max lattice speed<=.05, steps1..200000.");
    auto mesh=createMesh(c);Model::ParametersD p;Model model(p,mesh);prepareGeometry(model,c);prepareLattice(model,c);setInitialValues(model,c);simulate(model,c);return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
