// SPDX-License-Identifier: GPL-2.0-or-later
// Synthetic periodic FCC sphere bed: stock geometry, TRT, forcing and functors.
// Lifecycle follows installed laminar/poiseuille3d and porousMedia examples.
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>
using namespace olb;
using namespace olb::names;
using T=double;
using D=descriptors::D3Q27<descriptors::FORCE>;
using Model=Case<NavierStokes,Lattice<T,D>>;
using FlowDynamics=dynamics::Tuple<T,D,momenta::BulkTuple,equilibria::SecondOrder,collision::TRT,forcing::Kupershtokh>;
struct Config{int n=32;T gradient=10.,tau=.8,diameter=.0023;std::string warm;T warmScale=1.;};
constexpr T eps=.6,mu=.001,rho=1000.,magic=.1875;
T side(const Config& c){return c.diameter*std::cbrt(2*M_PI/(3*(1-eps)));}
Mesh<T,3> createMesh(const Config& c){
  T h=side(c)/c.n;IndicatorCuboid3D<T> box({side(c)-h,side(c)-h,side(c)-h},{h/2,h/2,h/2});
  Mesh<T,3> mesh(box,h,1);mesh.setOverlap(2);mesh.getCuboidDecomposition().setPeriodicity({true,true,true});return mesh;
}
void prepareGeometry(Model& m,const Config& c){
  auto& g=m.getGeometry();g.rename(0,1);const T l=side(c);
  // FCC centres and neighbouring images; all surfaces use the stock indicator.
  for(int i=-1;i<=3;++i)for(int j=-1;j<=3;++j)for(int k=-1;k<=3;++k){
    if((i+j+k)%2)continue;
    IndicatorSphere3D<T> sphere({i*l/2,j*l/2,k*l/2},c.diameter/2);g.rename(1,2,sphere);
  }
  g.communicate();g.checkForErrors();g.print();
  auto& b=g.getBlockGeometry(0);
  if(b.getNx()!=c.n||b.getNy()!=c.n||b.getNz()!=c.n)throw std::runtime_error("Unexpected periodic cell dimensions.");
}
void prepareLattice(Model& m,const Config& c){
  auto& l=m.getLattice(NavierStokes{});auto& g=m.getGeometry();
  l.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(c.n,c.tau,side(c),1e-4,mu/rho,rho);
  dynamics::set<FlowDynamics>(l,g,1);boundary::set<boundary::BounceBack>(l,g,2);
  l.setParameter<descriptors::OMEGA>(1/c.tau);l.setParameter<collision::TRT::MAGIC>(magic);
}
void setInitialValues(Model& m,const Config& c){
  auto& l=m.getLattice(NavierStokes{});auto& g=m.getGeometry();const auto& cv=l.getUnitConverter();
  T a=c.gradient/rho*cv.getPhysDeltaT()*cv.getPhysDeltaT()/cv.getPhysDeltaX();
  fields::set<descriptors::FORCE>(l,g.getMaterialIndicator({1}),Vector<T,3>{a,0.,0.});
  AnalyticalConst3D<T,T> rest(0.,0.,0.);momenta::setVelocity(l,g.getMaterialIndicator({1,2}),rest);
  if(!c.warm.empty()){
    std::ifstream stream(c.warm);std::string line;std::getline(stream,line);int count=0;
    auto& block=g.getBlockGeometry(0);
    while(std::getline(stream,line)){
      std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
      int i,j,k;T x,y,z,ux,uy,uz,density;row>>i>>j>>k>>x>>y>>z>>ux>>uy>>uz>>density;
      if(!row||i<0||j<0||k<0||i>=c.n||j>=c.n||k>=c.n||block.get({i,j,k})!=1)throw std::runtime_error("Invalid native initial field.");
      T velocity[3]={ux*c.warmScale*cv.getPhysDeltaT()/cv.getPhysDeltaX(),uy*c.warmScale*cv.getPhysDeltaT()/cv.getPhysDeltaX(),uz*c.warmScale*cv.getPhysDeltaT()/cv.getPhysDeltaX()};
      l.getBlock(0).get(i,j,k).iniEquilibrium(1.,velocity);++count;
    }
    if(count!=g.getStatistics().getNvoxel(1))throw std::runtime_error("Incomplete native initial field.");
  }
  l.initialize();
}
T meanVelocity(Model& m,const Config& c){
  auto& l=m.getLattice(NavierStokes{});auto& g=m.getGeometry();auto& b=g.getBlockGeometry(0);
  SuperLatticePhysVelocity3D<T,D> velocity(l,l.getUnitConverter());T sum=0;
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k)if(b.get({i,j,k})==1){
    int p[4]={0,i,j,k};T u[3];if(!velocity(u,p))throw std::runtime_error("Field read failed.");sum+=u[0];
  }
  return sum/(c.n*c.n*c.n);
}
void getResults(Model& m,const Config& c,int steps,T wall,bool converged){
  auto& l=m.getLattice(NavierStokes{});auto& g=m.getGeometry();auto& b=g.getBlockGeometry(0);const auto& cv=l.getUnitConverter();
  SuperLatticePhysVelocity3D<T,D> velocity(l,cv);SuperLatticeDensity3D<T,D> density(l);
  T maxu=0,maxdrho=0;int fluid=0;
  std::ofstream csv("field.csv");csv<<std::setprecision(17)<<"i,j,k,x_m,y_m,z_m,ux_m_s,uy_m_s,uz_m_s,density_kg_m3\n";
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k)if(b.get({i,j,k})==1){
    int p[4]={0,i,j,k};T u[3],r[1];if(!velocity(u,p)||!density(r,p))throw std::runtime_error("Field read failed.");
    auto x=b.getPhysR({i,j,k});csv<<i<<','<<j<<','<<k<<','<<x[0]<<','<<x[1]<<','<<x[2]<<','<<u[0]<<','<<u[1]<<','<<u[2]<<','<<r[0]*rho<<'\n';
    ++fluid;maxu=std::max(maxu,util::norm<3>(u));maxdrho=std::max(maxdrho,std::abs(r[0]-1));
  }
  std::ofstream raw("raw.json");raw<<std::setprecision(17)<<"{\"solver\":\"OpenLB 1.9r0\",\"collision\":\"TRT\",\"descriptor\":\"D3Q27\",\"forcing\":\"stock Kupershtokh\",\"magic\":"<<magic
    <<",\"resolution\":"<<c.n<<",\"particle_diameter_m\":"<<c.diameter<<",\"unit_cell_side_m\":"<<side(c)<<",\"declared_porosity\":"<<eps<<",\"fluid_nodes\":"<<fluid
    <<",\"gradient_Pa_m\":"<<c.gradient<<",\"viscosity_Pa_s\":"<<mu<<",\"density_kg_m3\":"<<rho<<",\"dt_s\":"<<cv.getPhysDeltaT()<<",\"dx_m\":"<<cv.getPhysDeltaX()
    <<",\"tau_even\":"<<c.tau<<",\"max_Mach\":"<<maxu*cv.getPhysDeltaT()/cv.getPhysDeltaX()*std::sqrt(3.)<<",\"max_density_departure\":"<<maxdrho
    <<",\"steps\":"<<steps<<",\"wall_s\":"<<wall<<",\"converged\":"<<(converged?"true":"false")<<"}\n";
}
void simulate(Model& m,const Config& c){
  auto& l=m.getLattice(NavierStokes{});auto start=std::chrono::steady_clock::now();
  util::ValueTracer<T> flow(6,1e-7,"superficial velocity"),energy(6,1e-7,"kinetic energy");
  int every=std::max(100,c.n*c.n/8);std::ofstream history("history.csv");history<<std::setprecision(17)<<"step,superficial_velocity_m_s,mean_lattice_energy\n";
  for(int i=1;i<=100000;++i){
    l.collideAndStream();if(i%every)continue;l.setProcessingContext(ProcessingContext::Evaluation);
    T u=meanVelocity(m,c),e=l.getStatistics().getAverageEnergy();
    if(!std::isfinite(u)||u*c.gradient<=0)throw std::runtime_error("Invalid native flow.");
    history<<i<<','<<u<<','<<e<<'\n';history.flush();flow.takeValue(u);energy.takeValue(e);
    const T observed=i*l.getUnitConverter().getPhysDeltaT()*(mu/rho)/(side(c)*side(c));
    if(flow.hasConverged()&&energy.hasConverged()&&observed>=.01){
      getResults(m,c,i,std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count(),true);return;
    }
    l.setProcessingContext(ProcessingContext::Simulation);
  }
  getResults(m,c,100000,std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count(),false);
  throw std::runtime_error("Native bed did not converge.");
}
int main(int argc,char** argv){
  initialize(&argc,&argv);try{
    if(argc!=5&&argc!=7)throw std::invalid_argument("Use: resolution gradient_Pa_m tau particle_diameter_m [native_initial_CSV scale]");
    Config c{std::stoi(argv[1]),std::stod(argv[2]),std::stod(argv[3]),std::stod(argv[4])};
    if(argc==7){c.warm=argv[5];c.warmScale=std::stod(argv[6]);}
    if(c.n<16||c.n>96||!std::isfinite(c.gradient)||c.gradient==0||!std::isfinite(c.tau)||c.tau<=.5||c.tau>1.5||!std::isfinite(c.diameter)||c.diameter<=0||!std::isfinite(c.warmScale))
      throw std::invalid_argument("Require N=16..96, nonzero finite gradient, tau (0.5,1.5], positive diameter.");
    auto mesh=createMesh(c);Model::ParametersD p;Model model(p,mesh);prepareGeometry(model,c);prepareLattice(model,c);setInitialValues(model,c);simulate(model,c);return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
