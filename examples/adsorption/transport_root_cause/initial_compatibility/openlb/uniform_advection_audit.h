// SPDX-License-Identifier: GPL-2.0-or-later
// Read-only one-step constant preservation diagnostic using stock equilibrium.
#ifndef P2P_UNIFORM_ADVECTION_AUDIT_H
#define P2P_UNIFORM_ADVECTION_AUDIT_H
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace p2p {
template<typename T, typename D, typename BLOCK, typename GEOMETRY>
void writeUniformAdvectionStencil(BLOCK& block, GEOMETRY& geometry,
                                 int nx,int ny,int nz,T dx,T dt,
                                 const std::string& path) {
  std::ofstream out(path);if(!out)throw std::runtime_error("Cannot write uniform advection stencil audit.");
  out<<std::setprecision(17);
  out<<"i,j,k,material,wall_adjacent_links,uniform_before_native,predicted_uniform_after_stream_native,uniform_change_native,scalar_equilibrium_divergence_per_s\n";
  auto wrap=[](int i,int n){return (i%n+n)%n;};
  for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k) {
    T next=T{1};int wallLinks=0;
    for(int p=0;p<D::q;++p) {
      const int cx=olb::descriptors::c<D>(p,0),cy=olb::descriptors::c<D>(p,1),cz=olb::descriptors::c<D>(p,2);
      const int ii=wrap(i-cx,nx),jj=wrap(j-cy,ny),kk=wrap(k-cz,nz);
      auto upstream=block.get(ii,jj,kk);
      auto velocity=upstream.template getField<olb::descriptors::VELOCITY>();
      next+=olb::equilibrium<D>::firstOrder(p,T{1},velocity.data());
      if((cx||cy||cz)&&geometry.get({ii,jj,kk})==2)++wallLinks;
    }
    // For initially uniform C=1, collision leaves first-order equilibrium
    // unchanged. Full-way BounceBack also leaves zero-u wall equilibrium
    // unchanged. The following is thus an exact first-step prediction.
    out<<i<<','<<j<<','<<k<<','<<geometry.get({i,j,k})<<','<<wallLinks<<",1,"<<next<<','<<next-T{1}<<','<<(T{1}-next)/dt<<'\n';
  }
}
}
#endif
