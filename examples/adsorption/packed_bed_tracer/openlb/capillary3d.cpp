// SPDX-License-Identifier: GPL-2.0-or-later
// Application diagnostics around unmodified OpenLB operators. Lifecycle follows
// advectionDiffusion1d and microMixer3d; FCC geometry follows our verified bed.
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>
using namespace olb;
using namespace olb::names;
using T=double;
using D=descriptors::D3Q27<descriptors::VELOCITY>;
using Model=Case<AdvectionDiffusion,Lattice<T,D>>;
struct Config { int n=24,cells=2; T side=.003993871093653167,Dm=1e-7,tauOdd=.7,end=150.;std::string field;bool uniform=false;T speed=.0005; };
constexpr T magic=.3;
int nx(const Config& c){return c.cells*c.n+1;}
T dx(const Config& c){return c.side/c.n;}
T dt(const Config& c){return (c.tauOdd-.5)*dx(c)*dx(c)/(descriptors::invCs2<T,D>()*c.Dm);}
Mesh<T,3> createMesh(const Config& c){
  const T h=dx(c);IndicatorCuboid3D<T> box({c.cells*c.side,c.side-h,c.side-h},{h/2,h/2,h/2});
  Mesh<T,3> mesh(box,h,1);mesh.setOverlap(2);mesh.getCuboidDecomposition().setPeriodicity({false,true,true});return mesh;
}
void prepareGeometry(Model& m,const Config& c){
  auto& g=m.getGeometry();T h=dx(c);
  IndicatorCuboid3D<T> segment({c.cells*c.side+h,c.side,c.side},{0.,0.,0.});g.rename(0,1,segment);
  if(!c.uniform) for(int i=-1;i<=2*c.cells+1;++i)for(int j=-1;j<=3;++j)for(int k=-1;k<=3;++k){
    if((i+j+k)%2)continue;
    // Same FCC sphere locations and radius as the accepted native flow field.
    T diameter=c.side/std::cbrt(2*M_PI/(3*.4));
    IndicatorSphere3D<T> sphere({i*c.side/2,j*c.side/2,k*c.side/2},diameter/2);g.rename(1,2,sphere);
  }
  IndicatorCuboid3D<T> inlet({h/2,c.side,c.side},{h/4,0.,0.});g.rename(1,3,inlet);
  IndicatorCuboid3D<T> outlet({h/2,c.side,c.side},{c.cells*c.side+h/4,0.,0.});g.rename(1,4,outlet);
  // In 1.9r0 indicator-based rename does not mark the super-geometry halo
  // dirty. The ordinary stock rename does so without changing these labels.
  g.rename(1,1);g.communicate();g.checkForErrors();g.print();
  auto& b=g.getBlockGeometry(0);
  if(b.getNx()!=nx(c)||b.getNy()!=c.n||b.getNz()!=c.n)throw std::runtime_error("Unexpected tracer mesh dimensions.");
}
void prepareLattice(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry();
  auto fi=g.getMaterialIndicator(1),oi=g.getMaterialIndicator(0);
  auto& bg=g.getBlockGeometry(0);
  for(int i:{0,nx(c)-1})for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k)if(bg.get({i,j,k})!=2){
    auto [kind,normal]=computeBoundaryTypeAndNormal(fi->getBlockIndicatorF(0),oi->getBlockIndicatorF(0),Vector<int,3>{i,j,k});
    if(util::norm<3>(normal)==0)throw std::runtime_error("Null stock boundary normal at "+std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+" neighbors "+std::to_string(bg.get({i+1,j,k}))+","+std::to_string(bg.get({i-1,j,k}))+","+std::to_string(bg.get({i,j+1,k}))+","+std::to_string(bg.get({i,j-1,k})));
  }
  l.setUnitConverter<AdeUnitConverter<T,D>>(dx(c),dt(c),c.cells*c.side,c.speed,c.Dm,1.);
  dynamics::set<AdvectionDiffusionTRTdynamics>(l,g.getMaterialIndicator({1,3,4}));
  boundary::set<boundary::BounceBack>(l,g,2);
  for(int material:{3,4}) boundary::set<boundary::AdvectionDiffusionDirichlet,AdvectionDiffusionTRTdynamics>(l,g,material);
  l.setParameter<descriptors::OMEGA>(1/(.5+magic/(c.tauOdd-.5)));
  l.setParameter<collision::TRT::MAGIC>(magic);
}
void setInitialValues(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry();auto& b=g.getBlockGeometry(0);
  std::vector<Vector<T,3>> u(c.n*c.n*c.n,Vector<T,3>{0.,0.,0.});std::vector<bool> seen(u.size(),false);
  if(!c.uniform){
    std::ifstream stream(c.field);if(!stream)throw std::runtime_error("Cannot read native velocity field.");
    std::string line;std::getline(stream,line);
    while(std::getline(stream,line)){
      std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
      int i,j,k;T x,y,z,ux,uy,uz,rho;row>>i>>j>>k>>x>>y>>z>>ux>>uy>>uz>>rho;
      if(!row||i<0||i>=c.n||j<0||j>=c.n||k<0||k>=c.n)throw std::runtime_error("Invalid native field index.");
      int index=(i*c.n+j)*c.n+k;if(seen[index])throw std::runtime_error("Duplicate native field index.");
      seen[index]=true;u[index]={ux,uy,uz};
    }
  }
  for(int i=0;i<nx(c);++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k){
    int index=((i%c.n)*c.n+j)*c.n+k,mat=b.get({i,j,k});
    if(!c.uniform && (mat!=2)!=seen[index])throw std::runtime_error("Native geometry/velocity mismatch.");
    auto velocity=c.uniform ? Vector<T,3>{c.speed,0.,0.}:u[index];velocity*=dt(c)/dx(c);
    auto cell=l.getBlock(0).get(i,j,k);cell.setField<descriptors::VELOCITY>(velocity);cell.iniEquilibrium(0.,velocity.data());
  }
  // Dirichlet momenta carry their prescribed density separately from populations.
  AnalyticalConst3D<T,T> clean(0.);momenta::setDensity(l,g.getMaterialIndicator({3,4}),clean);
  l.initialize();
  l.setStatisticsOff();
}
// Diagnostics of the stock streaming links crossing a fixed control surface.
// They are not a replacement flux scheme or a boundary implementation.
T linkFlux(Model& m,const Config& c,int left){
  auto& block=m.getLattice(AdvectionDiffusion{}).getBlock(0);T sum=0;
  for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k)for(int p=0;p<D::q;++p){
    if(descriptors::c<D>(p,0)==1){
      int jj=(j+descriptors::c<D>(p,1)+c.n)%c.n,kk=(k+descriptors::c<D>(p,2)+c.n)%c.n;
      sum+=block.get(left,j,k)[p]-block.get(left+1,jj,kk)[descriptors::opposite<D>(p)];
    }
  }
  return sum*std::pow(dx(c),3)/dt(c);
}
void getResults(Model& m,const Config& c,int step,T incoming,T outgoing,std::ofstream& history){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& b=m.getGeometry().getBlockGeometry(0);
  T fluid=0,wall=0,minc=1,maxc=0;
  for(int i=1;i<nx(c)-1;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k){
    auto cell=l.getBlock(0).get(i,j,k);T value=lbm<D>::computeRho(cell);
    if(b.get({i,j,k})==2)wall+=value;else{fluid+=value;minc=std::min(minc,value);maxc=std::max(maxc,value);}
  }
  T volume=std::pow(dx(c),3);fluid*=volume;wall*=volume;
  history<<step<<','<<step*dt(c)<<','<<fluid<<','<<wall<<','<<incoming<<','<<outgoing<<','<<fluid+wall-incoming+outgoing<<','<<minc<<','<<maxc<<'\n';
}
void simulate(Model& m,const Config& c){
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry();auto& b=g.getBlockGeometry(0);
  std::ofstream history("inventory.csv"),flux("flux.csv");history<<std::setprecision(17);flux<<std::setprecision(17);
  history<<"step,time_s,fluid_tracer_kg,wall_population_tracer_kg,cumulative_inlet_kg,cumulative_outlet_kg,mass_residual_kg,min_concentration_kg_m3,max_concentration_kg_m3\n";
  flux<<"time_s,inlet_flux_kg_s,outlet_flux_kg_s\n";
  getResults(m,c,0,0.,0.,history);T incoming=0,outgoing=0;
  AnalyticalConst3D<T,T> feed(1.);momenta::setDensity(l,g.getMaterialIndicator(3),feed);
  const int steps=std::ceil(c.end/dt(c));auto started=std::chrono::steady_clock::now();
  for(int s=1;s<=steps;++s){
    l.collide(); // stock public split, as needed to inspect actual postcollision links
    T ji=linkFlux(m,c,0),jo=linkFlux(m,c,nx(c)-2);
    incoming+=ji*dt(c);outgoing+=jo*dt(c);flux<<(s-.5)*dt(c)<<','<<ji<<','<<jo<<'\n';
    l.AndStream();
    if(s%10==0||s==steps)getResults(m,c,s,incoming,outgoing,history);
  }
  T fluidVolume=0,flow=0,maxu=0;
  std::ofstream field("final.csv");field<<std::setprecision(17)<<"i,j,k,material,x_m,y_m,z_m,concentration_kg_m3,ux_m_s,uy_m_s,uz_m_s\n";
  for(int i=0;i<nx(c);++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.n;++k){
    int mat=b.get({i,j,k});auto cell=l.getBlock(0).get(i,j,k);auto u=cell.getField<descriptors::VELOCITY>();
    auto x=b.getPhysR({i,j,k});T val=lbm<D>::computeRho(cell);
    if(mat!=2){
      maxu=std::max(maxu,util::norm<3>(u));
      if(i>0&&i<nx(c)-1){fluidVolume+=std::pow(dx(c),3);flow+=u[0]*dx(c)/dt(c)*std::pow(dx(c),2)/(nx(c)-2);}
    }
    field<<i<<','<<j<<','<<k<<','<<mat<<','<<x[0]<<','<<x[1]<<','<<x[2]<<','<<val<<','<<u[0]*dx(c)/dt(c)<<','<<u[1]*dx(c)/dt(c)<<','<<u[2]*dx(c)/dt(c)<<'\n';
  }
  std::ofstream raw("raw.json");raw<<std::setprecision(17)<<"{\"solver\":\"OpenLB 1.9r0\",\"collision\":\"stock AdvectionDiffusionTRTdynamics\",\"descriptor\":\"stock D3Q27<VELOCITY>\",\"boundary\":\"stock TRT Dirichlet step inlet / absorbing outlet\",\"n\":"<<c.n<<",\"cells\":"<<c.cells<<",\"nx\":"<<nx(c)<<",\"uniform\":"<<(c.uniform?"true":"false")<<",\"side_m\":"<<c.side<<",\"length_m\":"<<c.cells*c.side<<",\"area_m2\":"<<c.side*c.side<<",\"dx_m\":"<<dx(c)<<",\"dt_s\":"<<dt(c)<<",\"diffusivity_m2_s\":"<<c.Dm<<",\"tau_odd\":"<<c.tauOdd<<",\"tau_even\":"<<.5+magic/(c.tauOdd-.5)<<",\"magic\":"<<magic<<",\"steps\":"<<steps<<",\"fluid_volume_m3\":"<<fluidVolume<<",\"flow_m3_s\":"<<flow<<",\"max_lattice_speed\":"<<maxu<<",\"wall_s\":"<<std::chrono::duration<T>(std::chrono::steady_clock::now()-started).count()<<"}\n";
}
int main(int argc,char** argv){
  initialize(&argc,&argv);try{
    if(argc!=9)throw std::invalid_argument("Use: N cells side_m D_m2_s tau_odd end_s native_field_or_uniform uniform_speed_m_s");
    Config c;c.n=std::stoi(argv[1]);c.cells=std::stoi(argv[2]);c.side=std::stod(argv[3]);c.Dm=std::stod(argv[4]);c.tauOdd=std::stod(argv[5]);c.end=std::stod(argv[6]);c.field=argv[7];c.uniform=c.field=="uniform";c.speed=std::stod(argv[8]);
    if(c.n<16||c.n>64||c.cells<1||c.cells>8||!std::isfinite(c.side)||c.side<=0||!std::isfinite(c.Dm)||c.Dm<=0||!std::isfinite(c.tauOdd)||c.tauOdd<=.5||c.tauOdd>1.5||!std::isfinite(c.end)||c.end<=0||!std::isfinite(c.speed)||c.speed<=0||c.end/dt(c)>200000)
      throw std::invalid_argument("Require N=16..64, cells=1..8, positive finite SI scales, tau_odd in (.5,1.5], <=200000 steps.");
    auto mesh=createMesh(c);Model::ParametersD p;Model model(p,mesh);prepareGeometry(model,c);prepareLattice(model,c);setInitialValues(model,c);simulate(model,c);return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
