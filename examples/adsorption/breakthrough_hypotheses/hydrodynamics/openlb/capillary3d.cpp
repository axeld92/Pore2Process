// SPDX-License-Identifier: GPL-2.0-or-later
// Wall-bounded sphere packing with stock OpenLB geometry, TRT and forcing.
// Domain construction and simulation follow laminar/poiseuille3d and the
// project's external periodic-packing example.
#include <olb.h>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace olb;
using namespace olb::names;
using T = double;
using D = descriptors::D3Q27<descriptors::FORCE>;
using Model = Case<NavierStokes,Lattice<T,D>>;
using FlowDynamics = dynamics::Tuple<T,D,momenta::BulkTuple,equilibria::SecondOrder,
                                     collision::TRT,forcing::Kupershtokh>;
constexpr T diameter = .0008, viscosity = .001/1.2, density = 1000.;

struct Config {
  int nx, nt;
  T length, radius, tau, magic, throughput, targetInterstitial=0.;
  std::vector<Vector<T,3>> centres;
  std::string warmField;
  T warmDensityScale=1.;
  T spacing() const { return length/nx; }
  T centre() const { return nt*spacing()/2; }
};

Mesh<T,3> createMesh(const Config& c) {
  const T h = c.spacing();
  IndicatorCuboid3D<T> box({(c.nx-1)*h,(c.nt-1)*h,(c.nt-1)*h},{h/2,h/2,h/2});
  Mesh<T,3> mesh(box,h,1);
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,false,false});
  return mesh;
}

void prepareGeometry(Model& model, Config& c) {
  auto& g = model.getGeometry();
  const T h=c.spacing(), centre=c.centre();
  IndicatorCuboid3D<T> box({c.length,c.nt*h,c.nt*h},{0.,0.,0.});
  g.rename(0,3,box);
  IndicatorCylinder3D<T> tube({-h,centre,centre},{c.length+h,centre,centre},c.radius);
  g.rename(3,1,tube);
  for (const auto& centre : c.centres) {
    for (int i=-1; i<=1; ++i) {
      const Vector<T,3> image{centre[0]+i*c.length,centre[1],centre[2]};
      if (image[0]+diameter/2 >= 0. && image[0]-diameter/2 <= c.length) {
        IndicatorSphere3D<T> sphere(image,diameter/2);
        g.rename(1,2,sphere);
      }
    }
  }
  g.rename(1,1);
  g.communicate();
  g.checkForErrors();
  g.print();
  const auto& b = g.getBlockGeometry(0);
  if (b.getNx()!=c.nx || b.getNy()!=c.nt || b.getNz()!=c.nt) {
    throw std::runtime_error("Unexpected confined packing dimensions");
  }
  const T discreteFlowArea = g.getStatistics().getNvoxel(1)*h*h/c.nx;
  c.targetInterstitial=c.throughput/discreteFlowArea;
}

void prepareLattice(Model& model, const Config& c) {
  auto& l = model.getLattice(NavierStokes{});
  auto& g = model.getGeometry();
  l.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(
    c.nx,c.tau,c.length,c.targetInterstitial,viscosity/density,density);
  dynamics::set<FlowDynamics>(l,g,1);
  boundary::set<boundary::BounceBack>(l,g,2);
  boundary::set<boundary::BounceBack>(l,g,3);
  l.setParameter<descriptors::OMEGA>(1/c.tau);
  l.setParameter<collision::TRT::MAGIC>(c.magic);
}

void setDrivingGradient(Model& model, const Vector<T,3>& gradient) {
  auto& l = model.getLattice(NavierStokes{});
  auto& g = model.getGeometry();
  const auto& cv = l.getUnitConverter();
  const T f = cv.getPhysDeltaT()*cv.getPhysDeltaT()/(density*cv.getPhysDeltaX());
  fields::set<descriptors::FORCE>(l,g.getMaterialIndicator({1}),gradient*f);
}

void setInitialValues(Model& model, const Config& c, const Vector<T,3>& gradient) {
  auto& l = model.getLattice(NavierStokes{});
  auto& g = model.getGeometry();
  setDrivingGradient(model,gradient);
  AnalyticalConst3D<T,T> rest(0.,0.,0.);
  momenta::setVelocity(l,g.getMaterialIndicator({1,2,3}),rest);
  if (!c.warmField.empty()) {
    std::ifstream field(c.warmField,std::ios::binary);
    auto& block = g.getBlockGeometry(0);
    const auto& cv = l.getUnitConverter();
    const T velocityScale = cv.getPhysDeltaT()/cv.getPhysDeltaX();
    int rows=0;
    std::array<T,7> row;
    while (field.read(reinterpret_cast<char*>(row.data()),sizeof(row))) {
      const int i=int(row[0]),j=int(row[1]),k=int(row[2]);
      if (i<0 || i>=c.nx || j<0 || j>=c.nt || k<0 || k>=c.nt ||
          block.get({i,j,k})!=1) throw std::runtime_error("Invalid warm field index");
      const T rho = 1.+c.warmDensityScale*(row[6]/density-1.);
      T u[3]{row[3]*velocityScale,row[4]*velocityScale,row[5]*velocityScale};
      l.getBlock(0).get(i,j,k).iniEquilibrium(rho,u);
      ++rows;
    }
    if (rows!=g.getStatistics().getNvoxel(1) || field.gcount()!=0) {
      throw std::runtime_error("Incomplete warm field");
    }
  }
  l.initialize();
}

Vector<T,3> meanVelocity(Model& model, const Config& c) {
  auto& l = model.getLattice(NavierStokes{});
  auto& b = model.getGeometry().getBlockGeometry(0);
  SuperLatticePhysVelocity3D<T,D> velocity(l,l.getUnitConverter());
  Vector<T,3> total{0.,0.,0.};
  int fluid = 0;
  for (int i=0; i<c.nx; ++i) for (int j=0; j<c.nt; ++j) for (int k=0; k<c.nt; ++k) {
    if (b.get({i,j,k}) != 1) continue;
    int p[4]{0,i,j,k}; T u[3];
    if (!velocity(u,p)) throw std::runtime_error("Velocity read failed");
    for (unsigned a=0; a<3; ++a) total[a] += u[a];
    ++fluid;
  }
  return total*(1./fluid);
}

void getResults(Model& model, const Config& c, const Vector<T,3>& gradient, int steps, T wall) {
  auto& l = model.getLattice(NavierStokes{});
  auto& b = model.getGeometry().getBlockGeometry(0);
  const auto& cv = l.getUnitConverter();
  SuperLatticePhysVelocity3D<T,D> velocity(l,cv);
  SuperLatticeDensity3D<T,D> rho(l);
  std::ofstream field("field.f64",std::ios::binary);
  T maxSpeed=0., maxDensityError=0.;
  int fluid=0;
  for (int i=0; i<c.nx; ++i) for (int j=0; j<c.nt; ++j) for (int k=0; k<c.nt; ++k) {
    if (b.get({i,j,k}) != 1) continue;
    int p[4]{0,i,j,k}; T u[3],r[1];
    if (!velocity(u,p) || !rho(r,p)) throw std::runtime_error("Field read failed");
    const std::array<T,7> row{T(i),T(j),T(k),u[0],u[1],u[2],r[0]*density};
    field.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    maxSpeed=std::max(maxSpeed,util::norm<3>(u));
    maxDensityError=std::max(maxDensityError,std::abs(r[0]-1));
    ++fluid;
  }
  std::ofstream raw("raw.json");
  raw << std::setprecision(17)
      << "{\"solver\":\"OpenLB 1.9r0\",\"collision\":\"TRT\",\"descriptor\":\"D3Q27\","
      << "\"forcing\":\"stock Kupershtokh\",\"boundary\":\"stock BounceBack\",\"magic\":" << c.magic
      << ",\"shape\":[" << c.nx << ',' << c.nt << ',' << c.nt << ']'
      << ",\"length_m\":" << c.length << ",\"radius_m\":" << c.radius << ",\"tau_even\":" << c.tau
      << ",\"dx_m\":" << cv.getPhysDeltaX() << ",\"dt_s\":" << cv.getPhysDeltaT()
      << ",\"fluid_nodes\":" << fluid << ",\"gradient_Pa_m\":[" << gradient[0] << ',' << gradient[1] << ',' << gradient[2] << ']'
      << ",\"target_interstitial_m_s\":" << c.targetInterstitial
      << ",\"target_throughput_m3_s\":" << c.throughput
      << ",\"max_Mach\":" << maxSpeed*cv.getPhysDeltaT()/cv.getPhysDeltaX()*std::sqrt(3.)
      << ",\"max_density_departure\":" << maxDensityError << ",\"steps\":" << steps
      << ",\"wall_s\":" << wall << ",\"converged\":true}\n";
  if (!raw || !field) throw std::runtime_error("Cannot write flow field or diagnostics");
}

void simulate(Model& model, const Config& c, Vector<T,3> gradient) {
  auto& l = model.getLattice(NavierStokes{});
  const auto start = std::chrono::steady_clock::now();
  std::ofstream history("history.csv");
  history << std::setprecision(17) << "step,ux_m_s,uy_m_s,uz_m_s,lattice_energy,gx_Pa_m,gy_Pa_m,gz_Pa_m\n";
  const int every = 100;
  Vector<T,3> previous{0.,0.,0.};
  T previousEnergy=0.;
  int stable=0, adjusted=0;
  const T target = c.targetInterstitial;
  for (int step=1; step<=30000; ++step) {
    l.collideAndStream();
    if (step%every) continue;
    l.setProcessingContext(ProcessingContext::Evaluation);
    const auto u = meanVelocity(model,c);
    const T energy = l.getStatistics().getAverageEnergy();
    if (!std::isfinite(u[0]) || u[0]<=0. || !std::isfinite(energy) || energy<=0.) {
      throw std::runtime_error("Non-finite or reversed native flow");
    }
    history << step << ',' << u[0] << ',' << u[1] << ',' << u[2] << ',' << energy << ','
            << gradient[0] << ',' << gradient[1] << ',' << gradient[2] << '\n';
    history.flush();
    const T change = util::norm<3>(u-previous)/target;
    const T energyChange = std::abs(energy-previousEnergy)/energy;
    previous = u;
    previousEnergy = energy;
    const T residual = std::abs(target-u[0])/target;
    stable = change<1e-5 && energyChange<1e-5 ? stable+1 : 0;
    const T viscousTime = step*l.getUnitConverter().getPhysDeltaT()*(viscosity/density)/(diameter*diameter);
    if (stable>=3 && residual<=.002 && viscousTime>=.025) {
      const T elapsed = std::chrono::duration<T>(std::chrono::steady_clock::now()-start).count();
      getResults(model,c,gradient,step,elapsed);
      return;
    }
    l.setProcessingContext(ProcessingContext::Simulation);
    if (residual>.002 && change<3e-4 && energyChange<6e-4 && step-adjusted>=4*every) {
      gradient[0] += (target-u[0])*(.85*gradient[0]/u[0]);
      setDrivingGradient(model,gradient);
      stable=0;
      adjusted=step;
    }
  }
  throw std::runtime_error("Native flow failed the convergence or matched-throughput screen");
}

int main(int argc, char** argv) {
  initialize(&argc,&argv);
  try {
    if (argc!=10 && argc!=12) throw std::invalid_argument(
      "Use Nx Nt length_m radius_m tau magic throughput_m3_s gx_Pa_m centres [warm_field density_scale]");
    Config c{std::stoi(argv[1]),std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]),
             std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]),0.,{},argc==12 ? argv[10] : "",
             argc==12 ? std::stod(argv[11]) : 1.};
    const Vector<T,3> gradient{std::stod(argv[8]),0.,0.};
    if (c.nx<16 || c.nt<16 || !std::isfinite(c.length) || c.length<=diameter ||
        !std::isfinite(c.radius) || c.radius<=diameter || c.centre()<c.radius+c.spacing() ||
        !std::isfinite(c.tau) || c.tau<=.5 || c.tau>1.5 ||
        !std::isfinite(c.magic) || c.magic<=0. ||
        !std::isfinite(c.throughput) || c.throughput<=0. ||
        !std::isfinite(util::norm<3>(gradient)) || gradient[0]<=0. ||
        !std::isfinite(c.warmDensityScale) || c.warmDensityScale<=0.) {
      throw std::invalid_argument("Invalid native lattice or target flow");
    }
    std::ifstream centres(argv[9]);
    T x,y,z;
    while (centres >> x >> y >> z) c.centres.emplace_back(x,y,z);
    if (!centres.eof() || c.centres.empty()) {
      throw std::invalid_argument("Invalid or incomplete sphere-centre file");
    }
    auto mesh = createMesh(c);
    Model::ParametersD params;
    Model model(params,mesh);
    prepareGeometry(model,c);
    prepareLattice(model,c);
    setInitialValues(model,c,gradient);
    simulate(model,c,gradient);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
