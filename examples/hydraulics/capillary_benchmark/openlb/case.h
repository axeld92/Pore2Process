#pragma once
// Pore2Process hydraulic benchmark, using the OpenLB 1.9 API.
// Application case following examples/laminar/poiseuille3d/case.h.
// All dynamics, descriptors, geometry and integration functors are upstream.
// SPDX-License-Identifier: GPL-2.0-or-later
#include <olb.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>

using namespace olb;
using namespace olb::names;
using T = double;
using D = descriptors::D3Q19<>;
using Model = Case<NavierStokes, Lattice<T,D>>;

struct Config {
  int resolution = 16;
  int maxSteps = 40000;
  T diameter = 1e-4, length = 4e-4, cellWidth = 2e-4;
  T mu = 1e-3, rho = 1000., dp = 0.1, tau = 0.8, tolerance = 1e-8;
  bool geometryOnly = false;
};

struct Plane { T q=0, mass=0, pressure=0; int count=0; };


Mesh<T,3> createMesh(const Config& c) {
  const T dx=c.diameter/c.resolution;
  IndicatorCuboid3D<T> box({c.length+2*dx,c.diameter+dx,c.diameter+dx},
                           {-dx,-0.5*dx,-0.5*dx});
  Mesh<T,3> mesh(box,dx,1);
  mesh.setOverlap(3);
  return mesh;
}

void prepareGeometry(Model& model,const Config& c) {
  auto& geo=model.getGeometry();
  const T dx=c.diameter/c.resolution,radius=c.diameter/2;
  geo.rename(0,2);
  // Half-cell transverse centers avoid exact node/wall coincidences.
  // Use the stock cylinder and its stock distance calculation unchanged.
  IndicatorCylinder3D<T> fluid({-0.25*dx,radius,radius},
                               {c.length+0.25*dx,radius,radius}, radius);
  geo.rename(2,1,fluid);
  // Select whole lattice planes with half-cell-thick slabs, rather than placing
  // a zero-width or floating-point endpoint exactly on a geometry boundary.
  IndicatorCuboid3D<T> inlet({0.5*dx,c.diameter+2*dx,c.diameter+2*dx}, {-0.25*dx,-dx,-dx});
  IndicatorCuboid3D<T> outlet({0.5*dx,c.diameter+2*dx,c.diameter+2*dx}, {c.length-0.25*dx,-dx,-dx});
  geo.rename(1,3,inlet);
  geo.rename(1,4,outlet);
  geo.clean();
  geo.innerClean();
  geo.checkForErrors();
  geo.print();
  auto& stats = geo.getStatistics();
  int expected = 0;
  for (int j=0;j<c.resolution;++j) for (int k=0;k<c.resolution;++k) {
    const T y=(j+0.5)*dx-radius,z=(k+0.5)*dx-radius;
    if (y*y+z*z < std::pow(radius,2)) ++expected;
  }
  if (stats.getNvoxel(3) != expected || stats.getNvoxel(4) != expected ||
      std::abs(stats.getMinPhysR(3)[0]) > 0.01*dx ||
      std::abs(stats.getMaxPhysR(3)[0]) > 0.01*dx ||
      std::abs(stats.getMinPhysR(4)[0]-c.length) > 0.01*dx ||
      std::abs(stats.getMaxPhysR(4)[0]-c.length) > 0.01*dx)
    throw std::runtime_error("Inlet/outlet geometry count or physical location mismatch.");
}

void prepareLattice(Model& model,const Config& c) {
  auto& lattice=model.getLattice(NavierStokes{});
  auto& geo=model.getGeometry();
  const T dx=c.diameter/c.resolution,radius=c.diameter/2;
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(
      c.resolution,c.tau,c.diameter,1e-3,c.mu/c.rho,c.rho);
  const auto& converter=lattice.getUnitConverter();
  converter.print();
  dynamics::set<TRTdynamics>(lattice,geo,1);
  lattice.setParameter<collision::TRT::MAGIC>(T(3./16.));
  IndicatorCylinder3D<T> wall({-2*dx,radius,radius},{c.length+2*dx,radius,radius},radius);
  setBouzidiBoundary<T,D,BouzidiPostProcessor>(lattice,geo,2,wall,{1,3,4});
  boundary::set<boundary::LocalPressure>(lattice,geo,3);
  boundary::set<boundary::LocalPressure>(lattice,geo,4);
  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

void setInitialValues(Model& model,const Config& c) {
  auto& lattice=model.getLattice(NavierStokes{});
  auto& geo=model.getGeometry();
  // Zero velocity everywhere. Only the pressure initialization is linear;
  // no analytical velocity/flow/permeability enters the numerical solution.
  AnalyticalLinear3D<T,T> initialPressure(-c.dp/c.length,0,0,c.dp);
  AnalyticalConst3D<T,T> rest(0.,0.,0.);
  momenta::setPressureAndVelocity(lattice,geo.getMaterialIndicator({1,2,3,4}),initialPressure,rest);
  lattice.initialize();
}

// Application CSV diagnostics compose existing OpenLB field/plane functors.
// No custom integration, dynamics, descriptor or indicator implementation.
Plane getPlane(Model& model,const Config& c,int i,std::ostream* samples=nullptr) {
  auto& lattice=model.getLattice(NavierStokes{});
  auto& geo=model.getGeometry();
  const auto& converter=lattice.getUnitConverter();
  const T dx=c.diameter/c.resolution;
  SuperLatticePhysVelocity3D<T,D> velocity(lattice,converter);
  SuperLatticePhysPressure3D<T,D> pressure(lattice,converter);
  SuperLatticeDensity3D<T,D> density(lattice);
  auto hyperplane=Hyperplane3D<T>().originAt({i*dx,0.5*dx,0.5*dx})
                                  .spannedBy({0,1,0},{0,0,1});
  HyperplaneLattice3D<T> grid(hyperplane,dx,c.resolution,c.resolution);
  IndicatorCircle2D<T> circle({c.diameter/2-0.5*dx,c.diameter/2-0.5*dx},c.diameter/2);
  SuperPlaneIntegralFluxVelocity3D<T> volume(lattice,converter,geo,grid,
      geo.getMaterialIndicator({0,1,2,3,4}),
      circle,
      BlockDataReductionMode::Discrete);
  SuperPlaneIntegralFluxPressure3D<T> press(lattice,converter,geo,grid,
      geo.getMaterialIndicator({0,1,2,3,4}),
      circle,
      BlockDataReductionMode::Discrete);
  // Physical velocity times lattice density, scaled by the reference density.
  SuperPlaneIntegralFluxMass3D<T> mass(velocity,density,geo,c.rho,1.,grid,
      geo.getMaterialIndicator({0,1,2,3,4}),
      circle,
      BlockDataReductionMode::Discrete);
  T q[5]{},p[3]{},m[5]{};
  int input[1]{};
  if (!volume(q,input) || !press(p,input) || !mass(m,input) || q[1]<=0)
    throw std::runtime_error("Failed OpenLB plane integral.");
  Plane s{q[0],m[0],p[0]/p[1],int(std::llround(q[1]/(dx*dx)))};
  if (samples) {
    for (int j=0;j<c.resolution;++j) for (int k=0;k<c.resolution;++k) {
      const int point[4]={0,i+1,j+1,k+1};
      const int material=geo.get({0,i+1,j+1,k+1});
      if (material!=1 && material!=3 && material!=4) continue;
      T u[3],pr[1],rho[1];
      if (!velocity(u,point) || !pressure(pr,point) || !density(rho,point))
        throw std::runtime_error("Failed lattice field export.");
      *samples << i << ',' << j << ',' << k << ',' << i*dx << ',' << (j+0.5)*dx << ',' << (k+0.5)*dx
               << ',' << u[0] << ',' << u[1] << ',' << u[2] << ',' << pr[0] << ',' << rho[0]*c.rho << '\n';
    }
  }
  return s;
}

void getResults(Model& model,const Config& c,int step,int checkEvery,T qChange,T dpChange,
                std::chrono::steady_clock::time_point started) {
  const T dx=c.diameter/c.resolution;
  const int nx=int(std::llround(c.length/dx));
  const auto& converter=model.getLattice(NavierStokes{}).getUnitConverter();
  std::ofstream planes("planes.csv"),samples("midplane.csv");
  planes << std::setprecision(17) << "i,x_m,node_count,Q_m3_s,mass_flow_kg_s,mean_pressure_Pa\n";
  samples << std::setprecision(17) << "i,j,k,x_m,y_m,z_m,ux_m_s,uy_m_s,uz_m_s,pressure_Pa,density_kg_m3\n";
  for (int i=0;i<=nx;++i) {
    auto s=getPlane(model,c,i,i==nx/2 ? &samples : nullptr);
    planes << i << ',' << i*dx << ',' << s.count << ',' << s.q << ',' << s.mass << ',' << s.pressure << '\n';
  }
  const auto a=getPlane(model,c,0),b=getPlane(model,c,nx),m=getPlane(model,c,nx/2),p1=getPlane(model,c,nx/4),p2=getPlane(model,c,3*nx/4);
  T wallSeconds=std::chrono::duration<T>(std::chrono::steady_clock::now()-started).count();
  std::ofstream out("raw.json");
  out << std::setprecision(17)
      << "{\"schema_version\":2,\"solver\":\"OpenLB\",\"openlb_version\":\"1.9r0\",\"resolution\":" << c.resolution
      << ",\"axial_intervals\":" << nx << ",\"transverse_index_offset\":0.5,\"wall_method\":\"OpenLB IndicatorCylinder3D + BouzidiPostProcessor\""
      << ",\"diameter_m\":" << c.diameter << ",\"length_m\":" << c.length << ",\"cell_width_m\":" << c.cellWidth
      << ",\"viscosity_Pa_s\":" << c.mu << ",\"reference_density_kg_m3\":" << c.rho
      << ",\"imposed_delta_p_Pa\":" << c.dp << ",\"dx_m\":" << dx << ",\"dt_s\":" << converter.getPhysDeltaT()
      << ",\"tau\":" << c.tau << ",\"inlet_nodes\":" << a.count << ",\"outlet_nodes\":" << b.count
      << ",\"Q_in_m3_s\":" << a.q << ",\"Q_out_m3_s\":" << b.q << ",\"Q_mid_m3_s\":" << m.q
      << ",\"mass_in_kg_s\":" << a.mass << ",\"mass_out_kg_s\":" << b.mass
      << ",\"tap_length_m\":" << (3*nx/4-nx/4)*dx << ",\"tap_delta_p_Pa\":" << p1.pressure-p2.pressure
      << ",\"inlet_pressure_Pa\":" << a.pressure << ",\"outlet_pressure_Pa\":" << b.pressure
      << ",\"steps\":" << step << ",\"check_interval_steps\":" << checkEvery << ",\"steady_tolerance\":" << c.tolerance
      << ",\"convergence_method\":\"OpenLB ValueTracer relative standard deviation\",\"convergence_window_samples\":6"
      << ",\"final_relative_Q_change\":" << qChange << ",\"final_relative_dp_change\":" << dpChange
      << ",\"converged\":true,\"wall_s\":" << wallSeconds << "}\n";
}


void simulate(Model& model,const Config& c) {
  const auto started=std::chrono::steady_clock::now();
  const T dx=c.diameter/c.resolution;
  const int nx=int(std::llround(c.length/dx));
  auto& lattice=model.getLattice(NavierStokes{});
  const auto& converter=lattice.getUnitConverter();
  std::ofstream history("history.csv");
  history << std::setprecision(17) << "step,time_s,Q_in_m3_s,Q_mid_m3_s,Q_out_m3_s,mass_in_kg_s,mass_out_kg_s,delta_p_taps_Pa,relative_Q_change,relative_dp_change\n";
  const int checkEvery = std::max(100,c.resolution*c.resolution/2);
  util::ValueTracer<T> flowConvergence(6,c.tolerance,"midplane flow");
  util::ValueTracer<T> pressureConvergence(6,c.tolerance,"tap pressure drop");
  int step=0;
  T oldQ=0,oldDp=0,qChange=1,dpChange=1;
  bool converged=false;
  for (step=1;step<=c.maxSteps;++step) {
    lattice.collideAndStream();
    if (step%checkEvery) continue;
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    auto a=getPlane(model,c,0),m=getPlane(model,c,nx/2),b=getPlane(model,c,nx);
    T dp=getPlane(model,c,nx/4).pressure-getPlane(model,c,3*nx/4).pressure;
    qChange=std::abs(m.q-oldQ)/std::max(std::abs(m.q),1e-100);
    dpChange=std::abs(dp-oldDp)/std::max(std::abs(dp),1e-100);
    if (!std::isfinite(m.q) || !std::isfinite(dp)) throw std::runtime_error("Non-finite flow/pressure.");
    history << step << ',' << converter.getPhysTime(step) << ',' << a.q << ',' << m.q << ',' << b.q << ','
            << a.mass << ',' << b.mass << ',' << dp << ',' << qChange << ',' << dpChange << '\n';
    std::cout << "[hydraulic] step=" << step << " Q=" << m.q << " dP=" << dp
              << " Qchange=" << qChange << " dpChange=" << dpChange << std::endl;
    flowConvergence.takeValue(m.q);
    pressureConvergence.takeValue(dp);
    oldQ=m.q; oldDp=dp;
    if (flowConvergence.hasConverged() && pressureConvergence.hasConverged()) {
      converged=true;
      break;
    }
    lattice.setProcessingContext(ProcessingContext::Simulation);
  }
  lattice.setProcessingContext(ProcessingContext::Evaluation);
  if (!converged) throw std::runtime_error("Maximum steps reached before OpenLB flow/pressure ValueTracer convergence.");
  getResults(model,c,step,checkEvery,qChange,dpChange,started);
}
