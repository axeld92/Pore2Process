// SPDX-License-Identifier: GPL-2.0-or-later
// External case following OpenLB's advectionDiffusion1d and poiseuille3d examples.
// Physical inputs, application geometry and read-only diagnostics only.
#include <olb.h>
#include <bit>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
using namespace olb;
using namespace olb::names;
using T=double;
using D=descriptors::D3Q27<descriptors::VELOCITY>;
template<typename V,typename DD,typename M=momenta::AdvectionDiffusionBulkTuple>
using Transport=CombinedAdvectionDiffusionRLBdynamics<V,DD,AdvectionDiffusionTRTdynamics<V,DD>,M>;
using Model=Case<AdvectionDiffusion,Lattice<T,D>>;
constexpr T radius=.0004, alpha=3.0542369282271403, magic=.3;
struct Config {
  std::string kind;
  int n,nz,steps;
  T side,dt,diffusivity,speed,amplitude;
  T dx() const {return side/n;}
  T tauOdd() const {return .5+3*diffusivity*dt/(dx()*dx());}
  T tauEven() const {return .5+magic/(tauOdd()-.5);}
};

Mesh<T,3> createMesh(const Config& c) {
  const T h=c.dx();
  IndicatorCuboid3D<T> extent({c.side-h,c.side-h,(c.nz-1)*h},{h/2,h/2,h/2});
  Mesh<T,3> mesh(extent,h,1);mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,true,true});
  return mesh;
}

void prepareGeometry(Model& m,const Config& c) {
  auto& g=m.getGeometry();const T h=c.dx();
  IndicatorCuboid3D<T> extent({c.side,c.side,c.nz*h},{0.,0.,0.});
  g.rename(0,2,extent);
  IndicatorCylinder3D<T> cylinder({c.side/2,c.side/2,-h},{c.side/2,c.side/2,(c.nz+1)*h},radius);
  g.rename(2,1,cylinder);
  g.rename(1,1);g.communicate();g.checkForErrors();g.print();
  auto& b=g.getBlockGeometry(0);
  if(b.getNx()!=c.n||b.getNy()!=c.n||b.getNz()!=c.nz)throw std::runtime_error("Unexpected mesh dimensions");
}

void prepareLattice(Model& m,const Config& c) {
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry();
  l.setUnitConverter<AdeUnitConverter<T,D>>(c.dx(),c.dt,c.side,std::max(c.speed,T{1e-4}),c.diffusivity,1.);
  dynamics::set<Transport>(l,g,1);
  boundary::set<boundary::BounceBack>(l,g,2);
  l.setParameter<descriptors::OMEGA>(1/c.tauEven());
  l.setParameter<collision::TRT::MAGIC>(magic);
}

void setInitialValues(Model& m,const Config& c) {
  static_assert(std::endian::native==std::endian::little&&sizeof(T)==8&&std::numeric_limits<T>::is_iec559);
  auto& l=m.getLattice(AdvectionDiffusion{});auto& g=m.getGeometry().getBlockGeometry(0);
  const T omega0=3*std::sqrt(T{3})*c.speed/(2*radius);
  const T radialNorm=std::cyl_bessel_j(2.,alpha);
  std::ofstream mask("material.u8",std::ios::binary), velocity("velocity.f64",std::ios::binary);
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int iz=0;iz<c.nz;++iz) {
    auto x=g.getPhysR({i,j,iz});const int mat=g.get({i,j,iz});
    if(std::abs(x[0]-(i+.5)*c.dx())>1e-15||std::abs(x[1]-(j+.5)*c.dx())>1e-15||std::abs(x[2]-(iz+.5)*c.dx())>1e-15)
      throw std::runtime_error("Unexpected cell-centre coordinates");
    if(mat!=1&&mat!=2)throw std::runtime_error("Unexpected material");
    const T a=x[0]-c.side/2,b=x[1]-c.side/2,r=std::hypot(a,b);
    // Prescribed angular shear: peak speed U at r=R/sqrt(3), zero at r=R.
    const T omega=mat==1?omega0*std::max(T{},1-(r/radius)*(r/radius)):0.;
    Vector<T,3> u{-omega*b,omega*a,0.};
    const T value=.5+c.amplitude*std::cyl_bessel_j(2.,alpha*r/radius)/radialNorm*std::cos(2*std::atan2(b,a));
    mask.put(static_cast<char>(mat));
    u*=c.dt/c.dx();
    auto cell=l.getBlock(0).get(i,j,iz);cell.setField<descriptors::VELOCITY>(u);
    cell.iniEquilibrium(value,u.data());
  }
  l.initialize();l.setStatisticsOff();
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.nz;++k) {
    auto u=l.getBlock(0).get(i,j,k).getField<descriptors::VELOCITY>();
    const T physical[3]={u[0]*c.dx()/c.dt,u[1]*c.dx()/c.dt,u[2]*c.dx()/c.dt};
    velocity.write(reinterpret_cast<const char*>(physical),sizeof(physical));
  }
  if(!mask||!velocity)throw std::runtime_error("Cannot write initialized inputs");
}

// These are addresses for reading the active reflected populations. They do
// not participate in collision, streaming or boundary setup.
struct WallLink {int i,j,k,p;};
std::vector<WallLink> wallLinks(Model& m,const Config& c) {
  std::vector<WallLink> links;
  auto& g=m.getGeometry().getBlockGeometry(0);
  auto wrap=[](int i,int n){return (i+n)%n;};
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.nz;++k)if(g.get({i,j,k})==1) {
    for(int p=1;p<D::q;++p) {
      const int a=wrap(i+descriptors::c<D>(p,0),c.n),b=wrap(j+descriptors::c<D>(p,1),c.n),z=wrap(k+descriptors::c<D>(p,2),c.nz);
      if(g.get({a,b,z})==2)links.push_back({a,b,z,p});
    }
  }
  return links;
}

void capture(Model& m,const Config& c,std::vector<T>& values) {
  auto& b=m.getLattice(AdvectionDiffusion{}).getBlock(0);
  #pragma omp parallel for schedule(static)
  for(int i=0;i<c.n;++i)for(int j=0;j<c.n;++j)for(int k=0;k<c.nz;++k) {
    auto cell=b.get(i,j,k);values[(i*c.n+j)*c.nz+k]=lbm<D>::computeRho(cell);
  }
}

void writeField(const std::string& filename,const std::vector<T>& v) {
  std::ofstream f(filename,std::ios::binary);
  f.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(T));
  if(!f)throw std::runtime_error("Cannot write "+filename);
}

void simulate(Model& m,const Config& c) {
  auto& l=m.getLattice(AdvectionDiffusion{});auto& b=l.getBlock(0);auto& g=m.getGeometry().getBlockGeometry(0);
  const auto links=wallLinks(m,c);const T cellVolume=std::pow(c.dx(),3);
  std::ofstream stencil("stencil.csv");stencil<<std::setprecision(17)<<"cx,cy,cz,weight\n";
  for(int p=0;p<D::q;++p)stencil<<descriptors::c<D>(p,0)<<','<<descriptors::c<D>(p,1)<<','<<descriptors::c<D>(p,2)<<','<<descriptors::t<T,D>(p)<<'\n';
  std::vector<T> values(c.n*c.n*c.nz),peak(values.size());
  // The per-x reductions are accumulated in fixed order across thread counts.
  std::vector<T> mass(c.n),minima(c.n),maxima(c.n);
  std::ofstream h("inventory.csv");h<<std::setprecision(17)
    <<"step,time_s,fluid_inventory_kg,active_wall_population_kg,kinetic_inventory_kg,min_concentration_kg_m3,max_concentration_kg_m3\n";
  T largest=0;int peakStep=0;
  auto start=std::chrono::steady_clock::now();
  for(int step=0;step<=c.steps;++step) {
    if(step)l.collideAndStream();
    #pragma omp parallel for schedule(static)
    for(int i=0;i<c.n;++i) {
      T sum=0,lo=1e100,hi=-1e100;
      for(int j=0;j<c.n;++j)for(int k=0;k<c.nz;++k)if(g.get({i,j,k})==1) {
        auto cell=b.get(i,j,k);const T C=lbm<D>::computeRho(cell);
        if(!std::isfinite(C)){sum=std::numeric_limits<T>::quiet_NaN();continue;}
        sum+=C;lo=std::min(lo,C);hi=std::max(hi,C);
      }
      mass[i]=sum;minima[i]=lo;maxima[i]=hi;
    }
    T fluid=0,wall=0;for(T sum:mass)fluid+=sum;
    if(!std::isfinite(fluid))throw std::runtime_error("Non-finite fluid concentration");
    for(const auto& a:links) {
      auto cell=b.get(a.i,a.j,a.k);wall+=cell[a.p]+descriptors::t<T,D>(a.p);
    }
    if(!std::isfinite(wall))throw std::runtime_error("Non-finite active wall inventory");
    const T lo=*std::min_element(minima.begin(),minima.end()),hi=*std::max_element(maxima.begin(),maxima.end());
    h<<step<<','<<step*c.dt<<','<<fluid*cellVolume<<','<<wall*cellVolume<<','<<(fluid+wall)*cellVolume<<','<<lo<<','<<hi<<'\n';
    const T excursion=std::max(.5-c.amplitude-lo,hi-.5-c.amplitude);
    if(step==0||excursion>largest){largest=std::max(T{},excursion);peakStep=step;capture(m,c,peak);}
    if(step==0||step==c.steps/2||step==c.steps) {
      capture(m,c,values);writeField(step==0?"initial.f64":step==c.steps?"final.f64":"middle.f64",values);
    }
  }
  writeField("peak.f64",peak);
  std::ofstream f("raw.json");f<<std::setprecision(17)
    <<"{\"solver\":\"OpenLB 1.9r0\",\"collision\":\"stock CombinedAdvectionDiffusionRLBdynamics + AdvectionDiffusionTRTdynamics\",\"descriptor\":\"stock D3Q27 VELOCITY\",\"kind\":\""<<c.kind<<"\",\"n\":"<<c.n<<",\"nz\":"<<c.nz<<",\"steps\":"<<c.steps<<",\"side_m\":"<<c.side<<",\"dx_m\":"<<c.dx()<<",\"dt_s\":"<<c.dt<<",\"duration_s\":"<<c.steps*c.dt
    <<",\"radius_m\":"<<radius<<",\"speed_m_s\":"<<c.speed<<",\"diffusivity_m2_s\":"<<c.diffusivity<<",\"amplitude_kg_m3\":"<<c.amplitude<<",\"flow_law\":\"stationary circular wall, quadratic angular speed\""<<",\"tau_odd\":"<<c.tauOdd()<<",\"tau_even\":"<<c.tauEven()<<",\"magic\":"<<magic<<",\"active_wall_links\":"<<links.size()<<",\"peak_step\":"<<peakStep<<",\"maximum_excursion_kg_m3\":"<<largest<<",\"wall_s\":"<<std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count()<<"}\n";
  if(!f||!h)throw std::runtime_error("Cannot save diagnostics");
}

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    if(argc!=10)throw std::invalid_argument("Use kind N Nz steps L_m dt_s D_m2_s U_m_s amplitude_kg_m3");
    Config c{argv[1],std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]),std::stod(argv[8]),std::stod(argv[9])};
    if(c.kind!="stationary_disk"||c.n<16||c.n>128||c.nz<4||c.nz>c.n||c.steps<4||c.steps>50000)
      throw std::invalid_argument("Unsupported benchmark geometry or discretization");
    if(!std::isfinite(c.side)||c.side<=2*radius||!std::isfinite(c.dt)||c.dt<=0||!std::isfinite(c.diffusivity)||c.diffusivity<=0||!std::isfinite(c.speed)||c.speed<0||!std::isfinite(c.amplitude)||c.amplitude<0||c.amplitude>.25||c.speed*c.dt/c.dx()>.05)
      throw std::invalid_argument("Invalid SI inputs or lattice speed above .05");
    auto mesh=createMesh(c);Model::ParametersD p;Model m(p,mesh);
    prepareGeometry(m,c);prepareLattice(m,c);setInitialValues(m,c);simulate(m,c);return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;return 1;}
}
