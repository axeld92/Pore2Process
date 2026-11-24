// SPDX-License-Identifier: GPL-2.0-or-later
// Stock OpenLB application: rectangular, wall-bounded monolayer packing.
// Hydraulic x-periodic closure; the downstream conservative tracer has
// an open seam in a bead-free interval. No vendor code is changed.
#include <olb.h>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace olb;
using namespace olb::names;
using T = double;
using D = descriptors::D3Q19<descriptors::FORCE>;
using Model = Case<NavierStokes,Lattice<T,D>>;
using FlowDynamics = dynamics::Tuple<T,D,momenta::BulkTuple,equilibria::SecondOrder,
                                     collision::TRT,forcing::Kupershtokh>;
constexpr T magic = .1875;

struct Bead { Vector<T,3> centre; T radius; };
struct Config {
  int nx, ny, nz;
  T h, tau, flow, viscosity, density;
  std::vector<Bead> beads;
  std::string warmField;
  T warmTau=.58, warmGradient=0.;
};

Mesh<T,3> createMesh(const Config& c) {
  // Lateral solid cells at -h/2 and extent+h/2, so half-way
  // bounce-back walls enclose exactly ny*nz fluid-candidate cells.
  IndicatorCuboid3D<T> box({(c.nx-1)*c.h,(c.ny+1)*c.h,(c.nz+1)*c.h},
                           {c.h/2,-c.h/2,-c.h/2});
  Mesh<T,3> mesh(box,c.h,1);
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,false,false});
  return mesh;
}

void prepareGeometry(Model& model, const Config& c) {
  auto& g = model.getGeometry();
  g.rename(0,2);
  IndicatorCuboid3D<T> liquid({c.nx*c.h,c.ny*c.h,c.nz*c.h},{0.,0.,0.});
  g.rename(2,1,liquid);
  for (const auto& bead : c.beads) {
    IndicatorSphere3D<T> sphere(bead.centre,bead.radius);
    g.rename(1,2,sphere);
  }
  g.communicate();
  g.checkForErrors();
  g.print();
  const auto& b = g.getBlockGeometry(0);
  if (b.getNx()!=c.nx || b.getNy()!=c.ny+2 || b.getNz()!=c.nz+2) {
    throw std::runtime_error("Unexpected rectangular mesh dimensions");
  }
}

void prepareLattice(Model& model,const Config& c,T target) {
  auto& l = model.getLattice(NavierStokes{});
  auto& g = model.getGeometry();
  l.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(
    c.nz,c.tau,c.nz*c.h,target,c.viscosity/c.density,c.density);
  dynamics::set<FlowDynamics>(l,g,1);
  boundary::set<boundary::BounceBack>(l,g,2);
  l.setParameter<descriptors::OMEGA>(1/c.tau);
  l.setParameter<collision::TRT::MAGIC>(magic);
}

void setDrivingGradient(Model& model,const Config& c,T gradient) {
  auto& l = model.getLattice(NavierStokes{});
  auto& g = model.getGeometry();
  const auto& cv = l.getUnitConverter();
  const T f=cv.getPhysDeltaT()*cv.getPhysDeltaT()/(c.density*cv.getPhysDeltaX());
  fields::set<descriptors::FORCE>(l,g.getMaterialIndicator({1}),Vector<T,3>{gradient*f,0.,0.});
}

void setInitialValues(Model& model,const Config& c,T gradient) {
  auto& l=model.getLattice(NavierStokes{});
  auto& g=model.getGeometry();
  setDrivingGradient(model,c,gradient);
  AnalyticalConst3D<T,T> rest(0.,0.,0.);
  momenta::setVelocity(l,g.getMaterialIndicator({1,2}),rest);
  if(!c.warmField.empty()) {
    // A previously computed field is an initializer, never an accepted
    // new result. Scale the lattice density perturbation to the new dt.
    std::ifstream field(c.warmField,std::ios::binary);
    const auto& cv=l.getUnitConverter();
    const T velocityScale=cv.getPhysDeltaT()/cv.getPhysDeltaX();
    const T densityScale=std::pow((c.tau-.5)/(c.warmTau-.5),2);
    std::array<T,7> row;
    int read=0;
    auto& block=g.getBlockGeometry(0);
    while(field.read(reinterpret_cast<char*>(row.data()),sizeof(row))) {
      const int i=int(row[0]),j=int(row[1]),k=int(row[2]);
      if(i<0||i>=c.nx||j<1||j>c.ny||k<1||k>c.nz||block.get({i,j,k})!=1)
        throw std::runtime_error("Warm field does not match saved geometry");
      const T rho=1.+densityScale*(row[6]/c.density-1.);
      T u[3]{row[3]*velocityScale,row[4]*velocityScale,row[5]*velocityScale};
      l.getBlock(0).get(i,j,k).iniEquilibrium(rho,u);
      ++read;
    }
    if(read!=g.getStatistics().getNvoxel(1)||field.gcount()!=0)
      throw std::runtime_error("Incomplete warm field");
  }
  l.initialize();
}

struct FlowCheck {T mean=0., fieldChange=0., energy=0.;};
FlowCheck checkFlow(Model& model,const Config& c,T target,std::vector<T>& previous) {
  auto& l=model.getLattice(NavierStokes{});
  auto& b=model.getGeometry().getBlockGeometry(0);
  SuperLatticePhysVelocity3D<T,D> velocity(l,l.getUnitConverter());
  T total=0., delta=0., squared=0.;
  int fluid=0, sampled=0;
  for(int i=0;i<c.nx;++i) for(int j=1;j<=c.ny;++j) for(int k=1;k<=c.nz;++k) {
    if(b.get({i,j,k})!=1) continue;
    int p[4]{0,i,j,k}; T u[3];
    if(!velocity(u,p)) throw std::runtime_error("Velocity read failed");
    total+=u[0];
    if(fluid%101==0) {
      for(int a=0;a<3;++a) {
        const int index=3*sampled+a;
        if(index<int(previous.size())) {
          delta+=(u[a]-previous[index])*(u[a]-previous[index]);
          previous[index]=u[a];
        } else {
          delta+=u[a]*u[a];
          previous.push_back(u[a]);
        }
        squared+=u[a]*u[a];
      }
      ++sampled;
    }
    ++fluid;
  }
  return {total/fluid,std::sqrt(delta/std::max(squared,target*target)),
          l.getStatistics().getAverageEnergy()};
}

void getResults(Model& model,const Config& c,T gradient,int steps,T wall,T target) {
  auto& l=model.getLattice(NavierStokes{});
  auto& b=model.getGeometry().getBlockGeometry(0);
  const auto& cv=l.getUnitConverter();
  SuperLatticePhysVelocity3D<T,D> velocity(l,cv);
  SuperLatticeDensity3D<T,D> rho(l);
  std::ofstream field("field.f64",std::ios::binary);
  T maxSpeed=0., maxDensityError=0., total=0.;
  int fluid=0;
  for(int i=0;i<c.nx;++i) for(int j=1;j<=c.ny;++j) for(int k=1;k<=c.nz;++k) {
    if(b.get({i,j,k})!=1) continue;
    int p[4]{0,i,j,k}; T u[3],r[1];
    if(!velocity(u,p)||!rho(r,p)) throw std::runtime_error("Field read failed");
    const std::array<T,7> row{T(i),T(j),T(k),u[0],u[1],u[2],r[0]*c.density};
    field.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    total+=u[0];
    maxSpeed=std::max(maxSpeed,util::norm<3>(u));
    maxDensityError=std::max(maxDensityError,std::abs(r[0]-1.));
    ++fluid;
  }
  std::ofstream raw("raw.json");
  raw<<std::setprecision(17)
     <<"{\"solver\":\"OpenLB 1.9r0\",\"collision\":\"TRT\",\"descriptor\":\"D3Q19\","
     <<"\"forcing\":\"stock Kupershtokh\",\"magic\":"<<magic<<",\"shape_with_walls\":["
     <<c.nx<<','<<c.ny+2<<','<<c.nz+2<<"],\"dx_m\":"<<c.h
     <<",\"tau_even\":"<<c.tau<<",\"dt_s\":"<<cv.getPhysDeltaT()
     <<",\"fluid_nodes\":"<<fluid<<",\"gradient_Pa_m\":"<<gradient
     <<",\"flow_m3_s\":"<<total*c.h*c.h/c.nx
     <<",\"mean_interstitial_m_s\":"<<total/fluid
     <<",\"target_interstitial_m_s\":"<<target
     <<",\"max_Mach\":"<<maxSpeed*cv.getPhysDeltaT()/cv.getPhysDeltaX()*std::sqrt(3.)
     <<",\"max_density_departure\":"<<maxDensityError<<",\"steps\":"<<steps
     <<",\"wall_s\":"<<wall<<",\"converged\":true}\n";
}

void simulate(Model& model,const Config& c,T gradient,T target) {
  auto& l=model.getLattice(NavierStokes{});
  const auto started=std::chrono::steady_clock::now();
  std::ofstream history("history.csv");
  history<<std::setprecision(17)<<"step,mean_u_m_s,energy,gradient_Pa_m,field_change\n";
  const int every=std::max(100,c.nz*c.nz/2);
  T previous=0., previousEnergy=0.;
  int stable=0, adjusted=0;
  std::vector<T> previousField;
  for(int step=1;step<=40000;++step) {
    l.collideAndStream();
    if(step%every) continue;
    l.setProcessingContext(ProcessingContext::Evaluation);
    const auto f=checkFlow(model,c,target,previousField);
    if(!std::isfinite(f.mean)||f.mean<=0.||!std::isfinite(f.energy)||f.energy<=0.)
      throw std::runtime_error("Non-finite or reversed flow");
    history<<step<<','<<f.mean<<','<<f.energy<<','<<gradient<<','<<f.fieldChange<<'\n';
    history.flush();
    const T change=std::abs(f.mean-previous)/target;
    const T energyChange=std::abs(f.energy-previousEnergy)/f.energy;
    previous=f.mean; previousEnergy=f.energy;
    const T residual=std::abs(f.mean-target)/target;
    const T viscousTime=step*l.getUnitConverter().getPhysDeltaT()*(c.viscosity/c.density)/
                        std::pow(c.nz*c.h,2);
    stable=change<1e-5 && energyChange<2e-5 && f.fieldChange<2e-4 ? stable+1:0;
    if(stable>=3 && residual<=.002 && viscousTime>=.2) {
      getResults(model,c,gradient,step,std::chrono::duration<T>(
                 std::chrono::steady_clock::now()-started).count(),target);
      return;
    }
    l.setProcessingContext(ProcessingContext::Simulation);
    if(residual>.002 && change<5e-4 && energyChange<1e-3 && step-adjusted>=4*every) {
      gradient*=1.+.85*(target/f.mean-1.);
      setDrivingGradient(model,c,gradient);
      adjusted=step; stable=0;
    }
  }
  throw std::runtime_error("Native flow did not meet convergence and throughput criteria");
}

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    if(argc!=10 && argc!=13) throw std::invalid_argument(
      "Use nx ny nz h_m tau flow_m3_s mu_Pa_s rho_kg_m3 beads_m [warm_field old_tau old_gradient]");
    Config c{std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3]),std::stod(argv[4]),
             std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]),std::stod(argv[8]),{},
             argc==13 ? argv[10] : "",argc==13 ? std::stod(argv[11]) : .58,
             argc==13 ? std::stod(argv[12]) : 0.};
    if(c.nx<8||c.ny<8||c.nz<6||!std::isfinite(c.h)||c.h<=0.||c.tau<=.5||c.tau>1.2||
       !std::isfinite(c.flow)||c.flow<=0.||c.viscosity<=0.||c.density<=0.)
      throw std::invalid_argument("Invalid mesh or carrier inputs");
    std::ifstream centres(argv[9]);
    T x,y,z,r;
    while(centres>>x>>y>>z>>r) {
      if(r<=0.||x-r<0.||x+r>c.nx*c.h||y-r<-1e-12||y+r>c.ny*c.h+1e-12||
         z-r<-1e-12||z+r>c.nz*c.h+1e-12)
        throw std::invalid_argument("Bead outside rectangular extent");
      c.beads.push_back({{x,y,z},r});
    }
    if(!centres.eof()||c.beads.empty()) throw std::invalid_argument("Invalid bead file");
    auto mesh=createMesh(c);
    Model::ParametersD params;
    Model model(params,mesh);
    prepareGeometry(model,c);
    const T count=model.getGeometry().getStatistics().getNvoxel(1);
    const T target=c.flow*c.nx/(count*c.h*c.h);
    const T eps=count/(T(c.nx)*c.ny*c.nz);
    const T superficial=c.flow/(c.ny*c.nz*c.h*c.h);
    const T dp=2*c.beads.front().radius;
    const T gradient=c.warmField.empty() ?
      150*c.viscosity*superficial*(1-eps)*(1-eps)/(dp*dp*eps*eps*eps)
        +12*c.viscosity*superficial/std::pow(c.nz*c.h,2) : c.warmGradient;
    prepareLattice(model,c,target);
    setInitialValues(model,c,gradient);
    simulate(model,c,gradient,target);
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<e.what()<<std::endl;
    return 1;
  }
}
