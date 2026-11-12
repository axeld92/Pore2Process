// SPDX-License-Identifier: GPL-2.0-or-later
// Read-only diagnostics of stock OpenLB collision/streaming, not a flux scheme.
// Restricted to a single local block periodic in all three directions.
#ifndef P2P_POPULATION_AUDIT_H
#define P2P_POPULATION_AUDIT_H
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

namespace p2p {
template<typename T, typename DESCRIPTOR>
class PopulationAudit {
  struct Record {
    int material=0, wallLinks=0, negativePopulations=0;
    T before=0, post=0, predicted=0, toFluid=0, toWall=0;
    T minPopulation=0, momentum[3]={0,0,0};
  };
  int _nx,_ny,_nz;
  T _dx,_dt;
  std::vector<Record> _records;
  int index(int i,int j,int k) const { return (i*_ny+j)*_nz+k; }
  static int wrap(int i,int n) { return (i%n+n)%n; }
public:
  PopulationAudit(int nx,int ny,int nz,T dx,T dt)
    :_nx(nx),_ny(ny),_nz(nz),_dx(dx),_dt(dt),_records(nx*ny*nz) {
    if(nx<1||ny<1||nz<1||dx<=0||dt<=0)
      throw std::invalid_argument("Population audit requires positive mesh scales.");
  }
  template<typename BLOCK, typename GEOMETRY>
  void begin(BLOCK& block, GEOMETRY& geometry) {
    for(int i=0;i<_nx;++i)for(int j=0;j<_ny;++j)for(int k=0;k<_nz;++k) {
      auto cell=block.get(i,j,k); auto& r=_records[index(i,j,k)]; r=Record{};
      r.material=geometry.get({i,j,k});
      r.before=olb::lbm<DESCRIPTOR>::computeRho(cell);
    }
  }
  template<typename BLOCK, typename GEOMETRY>
  void afterCollision(BLOCK& block, GEOMETRY& geometry) {
    for(int i=0;i<_nx;++i)for(int j=0;j<_ny;++j)for(int k=0;k<_nz;++k) {
      auto cell=block.get(i,j,k);auto& r=_records[index(i,j,k)];
      r.post=olb::lbm<DESCRIPTOR>::computeRho(cell);r.predicted=T{1};
      r.minPopulation=std::numeric_limits<T>::infinity();
      for(int p=0;p<DESCRIPTOR::q;++p) {
        const int cx=olb::descriptors::c<DESCRIPTOR>(p,0);
        const int cy=olb::descriptors::c<DESCRIPTOR>(p,1);
        const int cz=olb::descriptors::c<DESCRIPTOR>(p,2);
        const int ii=wrap(i-cx,_nx),jj=wrap(j-cy,_ny),kk=wrap(k-cz,_nz);
        auto upstream=block.get(ii,jj,kk);
        r.predicted+=upstream[p];
        // OpenLB stores shifted populations; restore the quadrature weight
        // before testing positivity. The weight cancels in a signed link flux.
        const T population=cell[p]+olb::descriptors::t<T,DESCRIPTOR>(p);
        r.minPopulation=std::min(r.minPopulation,population);
        r.negativePopulations+=(population<0);
        r.momentum[0]+=cx*population;r.momentum[1]+=cy*population;r.momentum[2]+=cz*population;
        const int ni=wrap(i+cx,_nx),nj=wrap(j+cy,_ny),nk=wrap(k+cz,_nz);
        if(cx||cy||cz) {
          auto neighbour=block.get(ni,nj,nk);
          const T net=cell[p]-neighbour[olb::descriptors::opposite<DESCRIPTOR>(p)];
          if(geometry.get({ni,nj,nk})==2){r.toWall+=net;++r.wallLinks;}
          else r.toFluid+=net;
        }
      }
    }
  }
  template<typename BLOCK, typename GEOMETRY>
  void afterStream(BLOCK& block, GEOMETRY&, int step, const std::string& path) const {
    std::ofstream out(path);if(!out)throw std::runtime_error("Cannot write population audit.");
    out<<std::setprecision(17);
    out<<"step,i,j,k,material,wall_adjacent_links,concentration_before,concentration_after_collision,predicted_after_stream,concentration_after_stream,collision_delta,expected_stream_delta,local_balance_residual,net_to_fluid_per_s,net_to_wall_per_s,min_unshifted_post_population,negative_post_populations,post_momentum_x,post_momentum_y,post_momentum_z\n";
    const T volume=_dx*_dx*_dx;
    for(int i=0;i<_nx;++i)for(int j=0;j<_ny;++j)for(int k=0;k<_nz;++k) {
      auto cell=block.get(i,j,k);const auto& r=_records[index(i,j,k)];
      const T after=olb::lbm<DESCRIPTOR>::computeRho(cell);
      out<<step<<','<<i<<','<<j<<','<<k<<','<<r.material<<','<<r.wallLinks<<','
         <<r.before<<','<<r.post<<','<<r.predicted<<','<<after<<','<<r.post-r.before<<','
         <<r.predicted-r.post<<','<<after-r.predicted<<','<<r.toFluid*volume/_dt<<','
         <<r.toWall*volume/_dt<<','<<r.minPopulation<<','<<r.negativePopulations<<','
         <<r.momentum[0]<<','<<r.momentum[1]<<','<<r.momentum[2]<<'\n';
    }
  }
};
}
#endif
