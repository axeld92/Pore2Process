#pragma once
// SPDX-License-Identifier: GPL-2.0-or-later
// Project application following the FORCED branch of OpenLB 1.9r0
// examples/laminar/poiseuille3d/case.h and cylinder3d's stock STLreader.
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <stdexcept>
using namespace olb;
using namespace olb::names;
using T=double;
using D=descriptors::D3Q19<descriptors::FORCE>;
using Model=Case<NavierStokes,Lattice<T,D>>;
// Compose existing OpenLB facilities through its standard Tuple. This alias
// defines no new collision, forcing, momenta, descriptor or templated class.
// The stock Guo/ForcedTRT alias failed our circular oracle at unequal rates;
// the stock exact-difference forcing avoids that rate-specific combination.
using ChannelDynamics=dynamics::Tuple<T,D,momenta::BulkTuple,
    equilibria::SecondOrder,collision::TRT,forcing::Kupershtokh>;
constexpr const char* collisionName="TRT";
constexpr const char* forcingName="stock Kupershtokh exact difference";

struct Config {
  int resolution=24, axial=4, maxSteps=120000;
  T diameter=.0046, mu=.001, rho=1000, gradient=1., tau=.8, tolerance=1e-8,oddTau=1.6;
  bool circle=false;
  bool snapshots=false;
  std::string warm;
  T warmScale=1.;
};
T dx(const Config& c) { return c.diameter/c.resolution; }

std::shared_ptr<IndicatorF3D<T>> fluidShape(const Config& c) {
  const T h=c.diameter/2,d=dx(c);
  if (c.circle) return std::make_shared<IndicatorCylinder3D<T>>(
      Vector<T,3>{-10*d,h,h},Vector<T,3>{(c.axial+10)*d,h,h},h);
  return std::make_shared<STLreader<T>>("channel.stl",d,1.);
}

Mesh<T,3> createMesh(const Config& c) {
  const T d=dx(c);
  IndicatorCuboid3D<T> box({c.axial*d,c.diameter+d,c.diameter+d},{0.,-.5*d,-.5*d});
  Mesh<T,3> mesh(box,d,1);
  mesh.setOverlap(3);
  mesh.getCuboidDecomposition().setPeriodicity({true,false,false});
  return mesh;
}

void prepareGeometry(Model& model,const Config& c) {
  auto& geometry=model.getGeometry();
  geometry.rename(0,2);
  auto shape=fluidShape(c);
  if (c.circle) geometry.rename(2,1,*shape);
  else {
    // STLreader's default indication caches octree voxels. Use its existing
    // exact winding-number sign through the stock SDF indicator to classify
    // lattice centres against the declared triangles instead of a second mesh.
    auto stl=std::dynamic_pointer_cast<STLreader<T>>(shape);
    IndicatorSDF3D<T> exact([stl](Vector<T,3> p){return stl->template signedDistance<SignMode::EXACT>(p);});
    exact.getMin()=shape->getMin(); exact.getMax()=shape->getMax();
    geometry.rename(2,1,exact);
  }
  geometry.clean(); geometry.innerClean(); geometry.checkForErrors(); geometry.print();
  if (geometry.getStatistics().getNvoxel(1)==0) throw std::runtime_error("No channel fluid voxels.");
}

void prepareLattice(Model& model,const Config& c) {
  auto& lattice=model.getLattice(NavierStokes{}); auto& geometry=model.getGeometry();
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(
      c.resolution,c.tau,c.diameter,.001,c.mu/c.rho,c.rho);
  dynamics::set<ChannelDynamics>(lattice,geometry,1);
  // Same parameterization as stock multiComponent/bubbleChannel2d: a fixed
  // odd relaxation time, avoiding a divergent odd time as tau approaches .5.
  lattice.setParameter<collision::TRT::MAGIC>((c.tau-.5)*(c.oddTau-.5));
  // Stock halfway bounce-back is robust at the unresolved sharp STL tips.
  // Continuous geometry is approached by the explicitly checked voxel grids.
  boundary::set<boundary::BounceBack>(lattice,geometry,2);
  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
}

void setInitialValues(Model& model,const Config& c) {
  auto& lattice=model.getLattice(NavierStokes{}); auto& geometry=model.getGeometry();
  const auto& converter=lattice.getUnitConverter();
  const T force=c.gradient/c.rho*converter.getPhysDeltaT()*converter.getPhysDeltaT()/dx(c);
  fields::set<descriptors::FORCE>(lattice,geometry.getMaterialIndicator({1,2}),Vector<T,3>{force,0.,0.});
  AnalyticalConst3D<T,T> rest(0.,0.,0.);
  momenta::setVelocity(lattice,geometry.getMaterialIndicator({1,2}),rest);
  // Optional starting estimate from a prior native solution on this exact mesh.
  // Only stock cell equilibrium initialization is used. Every run subsequently
  // evolves and checks OpenLB convergence; no oracle is supplied as a solution.
  if (!c.warm.empty()) {
    std::ifstream stream(c.warm); std::string line; std::getline(stream,line);
    if (!stream) throw std::runtime_error("Cannot read native warm-start CSV.");
    auto& block=geometry.getBlockGeometry(0);
    int loaded=0;
    while (std::getline(stream,line)) {
      std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line);
      int j,k; T y,z,u; row>>j>>k>>y>>z>>u;
      if (!row) throw std::runtime_error("Invalid native warm-start row.");
      T velocity[3]={u*c.warmScale*converter.getPhysDeltaT()/dx(c),0.,0.};
      for (int i=0;i<block.getNx();++i) if (block.get({i,j,k})==1)
        lattice.getBlock(0).get(i,j,k).iniEquilibrium(1.,velocity);
      ++loaded;
    }
    if (!loaded) throw std::runtime_error("Empty warm start.");
  }
  lattice.initialize();
}

struct Plane {T q=0,mass=0,maxSpeed=0,maxDensity=0;};
Plane plane(Model& model,const Config& c,T x) {
  auto& lattice=model.getLattice(NavierStokes{}); auto& geometry=model.getGeometry();
  const auto& converter=lattice.getUnitConverter(); const T d=dx(c);
  SuperLatticePhysVelocity3D<T,D> velocity(lattice,converter);
  SuperLatticeDensity3D<T,D> density(lattice);
  HyperplaneLattice3D<T> grid(Hyperplane3D<T>().originAt({x,.5*d,.5*d}).spannedBy({0,1,0},{0,0,1}),d,c.resolution,c.resolution);
  IndicatorCuboid2D<T> square({c.diameter,c.diameter},{-.25*d,-.25*d});
  SuperPlaneIntegralFluxVelocity3D<T> volume(lattice,converter,geometry,grid,
      geometry.getMaterialIndicator({1}),square,BlockDataReductionMode::Discrete);
  SuperPlaneIntegralFluxMass3D<T> mass(velocity,density,geometry,c.rho,1.,grid,
      geometry.getMaterialIndicator({1}),square,BlockDataReductionMode::Discrete);
  T q[5]{},m[5]{}; int input[1]{};
  if (!volume(q,input)||!mass(m,input)||q[1]<=0) throw std::runtime_error("Failed stock plane integral.");
  return {q[0],m[0],0.,0.};
}

void getResults(Model& model,const Config& c,int step,int every,T change,T elapsed,bool converged=true) {
  auto& lattice=model.getLattice(NavierStokes{}); auto& geometry=model.getGeometry();
  auto& block=geometry.getBlockGeometry(0); const auto& converter=lattice.getUnitConverter();
  SuperLatticePhysVelocity3D<T,D> velocity(lattice,converter);
  SuperLatticeDensity3D<T,D> density(lattice);
  std::ofstream samples("section.csv");
  samples<<std::setprecision(17)<<"j,k,y_m,z_m,ux_m_s,uy_m_s,uz_m_s,density_kg_m3\n";
  T maxSpeed=0,maxDensity=0,transverse=0; int count=0;
  for (int j=0;j<block.getNy();++j) for (int k=0;k<block.getNz();++k) {
    if (block.get({1,j,k})!=1) continue;
    const int point[4]={0,1,j,k}; auto pos=block.getPhysR({1,j,k}); T u[3],rho[1];
    if (!velocity(u,point)||!density(rho,point)) throw std::runtime_error("Native field query failed.");
    samples<<j<<','<<k<<','<<pos[1]<<','<<pos[2]<<','<<u[0]<<','<<u[1]<<','<<u[2]<<','<<rho[0]*c.rho<<'\n';
    maxSpeed=std::max(maxSpeed,std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]));
    transverse=std::max(transverse,std::sqrt(u[1]*u[1]+u[2]*u[2]));
    maxDensity=std::max(maxDensity,std::abs(rho[0]-1)); ++count;
  }
  std::ofstream planes("planes.csv"); planes<<std::setprecision(17)<<"x_m,Q_m3_s,mass_kg_s\n";
  for (int i=1;i<c.axial;++i) {auto p=plane(model,c,i*dx(c)); planes<<i*dx(c)<<','<<p.q<<','<<p.mass<<'\n';}
  const auto p=plane(model,c,dx(c));
  std::ofstream raw("raw.json"); raw<<std::setprecision(17)
    <<"{\"solver\":\"OpenLB\",\"version\":\"1.9r0\",\"collision_model\":\""<<collisionName<<"\",\"forcing\":\""<<forcingName<<"\",\"shape\":\""<<(c.circle?"circle":"six_point_star")<<"\",\"resolution\":"<<c.resolution
    <<",\"periodic_cells_requested\":"<<c.axial<<",\"actual_axial_planes\":"<<block.getNx()
    <<",\"diameter_m\":"<<c.diameter<<",\"dx_m\":"<<dx(c)<<",\"dt_s\":"<<converter.getPhysDeltaT()
    <<",\"viscosity_Pa_s\":"<<c.mu<<",\"density_kg_m3\":"<<c.rho<<",\"gradient_Pa_m\":"<<c.gradient<<",\"tau\":"<<c.tau
    <<",\"odd_tau\":"<<c.oddTau
    <<",\"wall_method\":\"stock halfway BounceBack\""
    <<",\"Q_m3_s\":"<<p.q<<",\"section_nodes\":"<<count<<",\"max_speed_m_s\":"<<maxSpeed
    <<",\"max_transverse_speed_m_s\":"<<transverse<<",\"max_density_departure\":"<<maxDensity
    <<",\"max_Mach\":"<<maxSpeed*converter.getPhysDeltaT()/dx(c)*std::sqrt(3.)
    <<",\"steps\":"<<step<<",\"check_interval\":"<<every<<",\"final_relative_change\":"<<change
    <<",\"tolerance\":"<<c.tolerance<<",\"converged\":"<<(converged?"true":"false")<<",\"wall_s\":"<<elapsed<<"}\n";
}

void simulate(Model& model,const Config& c) {
  const auto start=std::chrono::steady_clock::now(); auto& lattice=model.getLattice(NavierStokes{});
  const int every=std::max(100,c.resolution*c.resolution/2);
  util::ValueTracer<T> tracer(6,c.tolerance,"native axial flow");
  std::ofstream history("history.csv"); history<<std::setprecision(17)<<"step,Q_m3_s,relative_change\n";
  T old=0,change=1;
  for (int step=1;step<=c.maxSteps;++step) {
    lattice.collideAndStream(); if (step%every) continue;
    lattice.setProcessingContext(ProcessingContext::Evaluation); auto p=plane(model,c,dx(c));
    change=std::abs(p.q-old)/std::max(std::abs(p.q),1e-100); old=p.q;
    if (!std::isfinite(p.q)||p.q*c.gradient<=0) throw std::runtime_error("Nonfinite or wrong-signed native flow.");
    history<<step<<','<<p.q<<','<<change<<'\n'; history.flush(); tracer.takeValue(p.q);
    if (tracer.hasConverged()) {
      getResults(model,c,step,every,change,std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count()); return;
    }
    if(c.snapshots && step%(every*10)==0)
      getResults(model,c,step,every,change,std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count(),false);
    lattice.setProcessingContext(ProcessingContext::Simulation);
  }
  throw std::runtime_error("Native ValueTracer did not converge before maximum steps.");
}
