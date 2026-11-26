// SPDX-License-Identifier: GPL-2.0-or-later
// Periodic water flow through the *saved* tube voxel mask. Julia passes the
// local feed fraction and receives water velocity. Concentration is dimensionless.
#include <olb.h>
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

using namespace olb;
using namespace olb::names;
using T = double;
using D = descriptors::D3Q27<descriptors::FORCE>;
using Model = Case<NavierStokes,Lattice<T,D>>;
using FlowDynamics = dynamics::Tuple<T,D,momenta::BulkTuple,equilibria::SecondOrder,
                                 collision::TRT,forcing::Kupershtokh>;
using Row = std::array<T,7>;
struct Config {
  int period, transverse, repeats;
  T h, nu, rho, speedScale, axialAcceleration, tau;
  int length() const { return period*repeats; }
};

auto createMesh(const Config& c) {
  IndicatorCuboid3D<T> box({(c.length()-1)*c.h,(c.transverse-1)*c.h,(c.transverse-1)*c.h},
                           {c.h/2,c.h/2,c.h/2});
  Mesh<T,3> mesh(box,c.h,1);
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,false,false});
  return mesh;
}
void prepareGeometry(Model& model, const Config& c, const std::vector<Row>& rows) {
  auto& g=model.getGeometry();
  g.setWriteIncrementalVTK(false);
  IndicatorCuboid3D<T> box({c.length()*c.h,c.transverse*c.h,c.transverse*c.h},{0.,0.,0.});
  g.rename(0,2,box);
  auto& b=g.getBlockGeometry(0);
  if (b.getNx()!=c.length() || b.getNy()!=c.transverse || b.getNz()!=c.transverse)
    throw std::runtime_error("Unexpected voxel grid");
  for (int r=0;r<c.repeats;++r) for (const Row& v : rows)
    b.set({int(v[0])+r*c.period,int(v[1]),int(v[2])},1);
  g.rename(1,1);
  g.communicate();
  if (g.getStatistics().getNvoxel(1)!=rows.size()*c.repeats)
    throw std::runtime_error("Saved and reconstructed pore volumes differ");
}
void prepareLattice(Model& model, const Config& c) {
  auto& l=model.getLattice(NavierStokes{});
  auto& g=model.getGeometry();
  l.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(
    c.period,c.tau,c.period*c.h,2.24e-5,c.nu,c.rho);
  dynamics::set<FlowDynamics>(l,g,1);
  boundary::set<boundary::BounceBack>(l,g,2);
  l.setParameter<descriptors::OMEGA>(1/c.tau);
  l.setParameter<collision::TRT::MAGIC>(3./16.);
}
void setInitialValues(Model& model, const Config& c, const std::vector<Row>& rows) {
  auto& l=model.getLattice(NavierStokes{});
  auto& g=model.getGeometry();
  const T dt=l.getUnitConverter().getPhysDeltaT(), scale=dt/c.h;
  fields::set<descriptors::FORCE>(l,g.getMaterialIndicator({1}),
                                 Vector<T,3>{c.axialAcceleration*dt*scale,0.,0.});
  AnalyticalConst3D<T,T> rest(0.,0.,0.);
  momenta::setVelocity(l,g.getMaterialIndicator({1,2}),rest);
  for (int r=0;r<c.repeats;++r) for (const Row& v : rows) {
    T u[3]{v[3]*scale/c.speedScale,v[4]*scale/c.speedScale,v[5]*scale/c.speedScale};
    l.getBlock(0).get(int(v[0])+r*c.period,int(v[1]),int(v[2])).iniEquilibrium(v[6]/1000.,u);
  }
  l.initialize();
}
void simulate(Model& model, const Config& c, const std::vector<Row>& rows) {
  auto& l=model.getLattice(NavierStokes{});
  auto& block=l.getBlock(0);
  const T dt=l.getUnitConverter().getPhysDeltaT(), factor=dt*dt/c.h;
  const size_t n=rows.size()*c.repeats;
  std::vector<T> salt(n);
  std::vector<std::array<T,3>> u(n);
  std::vector<size_t> ids;
  ids.reserve(n);
  for (int r=0;r<c.repeats;++r) for (const Row& v : rows)
    ids.push_back(block.getCellId({int(v[0])+r*c.period,int(v[1]),int(v[2])}));
  std::cout << "SALT_READY " << dt << std::endl;
  int steps;
  T gravity, pump;
  while (std::cin >> steps >> gravity >> pump) {
    if (steps<0) break;
    std::ifstream in("salt.bin",std::ios::binary);
    if (!in.read(reinterpret_cast<char*>(salt.data()),n*sizeof(T))) throw std::runtime_error("Incomplete salt field");
    #pragma omp parallel for
    for (size_t p=0;p<n;++p) {
      block.get(ids[p]).template setField<descriptors::FORCE>(Vector<T,3>{
        pump*c.axialAcceleration*factor,0.,gravity*salt[p]*factor});
    }
    l.setProcessingContext(ProcessingContext::Simulation);
    for (int iter=0;iter<steps;++iter) l.collideAndStream();
    l.setProcessingContext(ProcessingContext::Evaluation);
    T sum=0., maxSpeed2=0., maxRho=0.;
    #pragma omp parallel for reduction(+:sum) reduction(max:maxSpeed2,maxRho)
    for (size_t p=0;p<n;++p) {
      T rho;
      block.get(ids[p]).computeRhoU(rho,u[p].data());
      maxRho=std::max(maxRho,std::abs(rho-1));
      maxSpeed2=std::max(maxSpeed2,util::normSqr<T,3>(u[p].data()));
      for (T& part : u[p]) part*=c.h/dt;
      sum+=u[p][0];
    }
    std::ofstream out("water.bin",std::ios::binary | std::ios::trunc);
    if (!out.write(reinterpret_cast<const char*>(u.data()),n*sizeof(u[0]))) throw std::runtime_error("Water output failed");
    out.close();
    std::cout << std::setprecision(14) << "SALT_DONE " << sum*c.h*c.h/c.length()
              << ' ' << std::sqrt(3*maxSpeed2) << ' ' << maxRho << std::endl;
  }
}
int main(int argc,char** argv) {
  initialize(&argc,&argv);
  try {
    if (argc!=11) throw std::invalid_argument("Use nx nt repeats h nu rho nativeSpeedScale axialAcceleration tau field");
    Config c{std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3]),std::stod(argv[4]),std::stod(argv[5]),
             std::stod(argv[6]),std::stod(argv[7]),std::stod(argv[8]),std::stod(argv[9])};
    std::ifstream f(argv[10],std::ios::binary);
    std::vector<Row> rows;
    Row v;
    while (f.read(reinterpret_cast<char*>(v.data()),sizeof(v))) rows.push_back(v);
    if (rows.empty() || !f.eof() || f.gcount()!=0) throw std::runtime_error("Invalid saved pore field");
    auto mesh=createMesh(c);
    Model::ParametersD params;
    Model model(params,mesh);
    prepareGeometry(model,c,rows);
    prepareLattice(model,c);
    setInitialValues(model,c,rows);
    simulate(model,c,rows);
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << std::endl; return 1; }
}
