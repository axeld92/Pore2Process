// SPDX-License-Identifier: GPL-2.0-or-later
// Application configuration/read-only verification of stock Bouzidi operators.
#ifndef P2P_HALFWAY_SETUP_H
#define P2P_HALFWAY_SETUP_H
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace p2p {
template<typename T,typename D,typename LATTICE,typename GEOMETRY,typename INDICATOR>
void configureHalfway(LATTICE& lattice,GEOMETRY& geometry,INDICATOR& indicator,
                      int nx,int ny,int nz,bool periodicY) {
  olb::dynamics::set<olb::NoDynamics>(lattice,geometry,2);
  olb::setBouzidiBoundary<T,D,olb::BouzidiPostProcessor>(lattice,geometry,2,indicator);
  auto& block=lattice.getBlock(0);auto& g=geometry.getBlockGeometry(0);
  // Only configure links that the stock setter registered; never construct
  // a replacement postprocessor or change an inactive population direction.
  for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k)if(g.get({i,j,k})==1){
    auto cell=block.get(i,j,k);auto q=cell.template getFieldPointer<olb::descriptors::BOUZIDI_DISTANCE>();
    for(int p=1;p<D::q;++p)if(q[p]>=0)q[p]=T{.5};
  }
  std::size_t links=0,seams=0,unexpected=0,missing=0,inactiveChanged=0;
  std::ofstream stencil("stencil.csv"),mask("mask.csv"),fields("halfway_link_fields.csv");
  stencil<<std::setprecision(17)<<"p,opposite,cx,cy,cz,weight\n";
  for(int p=0;p<D::q;++p)stencil<<p<<','<<olb::descriptors::opposite<D>(p)<<','<<olb::descriptors::c<D>(p,0)<<','<<olb::descriptors::c<D>(p,1)<<','<<olb::descriptors::c<D>(p,2)<<','<<olb::descriptors::t<T,D>(p)<<'\n';
  mask<<"i,j,k,material\n";fields<<"i,j,k,p,q\n";
  for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k)mask<<i<<','<<j<<','<<k<<','<<g.get({i,j,k})<<'\n';
  auto wrap=[](int i,int n){return (i%n+n)%n;};
  for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k)if(g.get({i,j,k})==1){
    auto cell=block.get(i,j,k);auto q=cell.template getFieldPointer<olb::descriptors::BOUZIDI_DISTANCE>();
    for(int p=1;p<D::q;++p){
      if(q[p]>=0)fields<<i<<','<<j<<','<<k<<','<<p<<','<<q[p]<<'\n';
      int ii=i+olb::descriptors::c<D>(p,0),jj=j+olb::descriptors::c<D>(p,1),kk=k+olb::descriptors::c<D>(p,2);
      bool seam=ii<0||ii>=nx||kk<0||kk>=nz||(periodicY&&(jj<0||jj>=ny));
      int ni=wrap(ii,nx),nj=periodicY?wrap(jj,ny):jj,nk=wrap(kk,nz);
      bool target=nj>=0&&nj<ny&&g.get({ni,nj,nk})==2;
      if(target){++links;if(seam)++seams;if(q[p]!=T{.5})++missing;}
      else {if(q[p]>=0)++unexpected;if(q[p]!=T{-1})++inactiveChanged;}
    }
  }
  std::ofstream out("halfway_links.json");out<<"{\"fluid_solid_links\":"<<links<<",\"periodic_seam_links\":"<<seams<<",\"missing_links\":"<<missing<<",\"unexpected_active_links\":"<<unexpected<<",\"inactive_changed\":"<<inactiveChanged<<",\"q\":0.5,\"operator\":\"stock BouzidiPostProcessor at PostStream\"}\n";
  if(!links||missing||unexpected||inactiveChanged)throw std::runtime_error("Stock half-way link coverage failed.");
}

template<typename T,typename D>
class HalfwayAudit {
  int nx,ny,nz;bool py;std::vector<T> predicted;
  int index(int i,int j,int k)const{return(i*ny+j)*nz+k;}
  static int wrap(int i,int n){return(i%n+n)%n;}
public:
  HalfwayAudit(int x,int y,int z,bool periodicY):nx(x),ny(y),nz(z),py(periodicY),predicted(x*y*z){}
  template<typename BLOCK,typename GEO>
  void beforeStream(BLOCK& b,GEO& g){
    for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k)if(g.get({i,j,k})==1){
      T value=1;auto cell=b.get(i,j,k);
      for(int p=0;p<D::q;++p){
        int ii=wrap(i-olb::descriptors::c<D>(p,0),nx),jj=j-olb::descriptors::c<D>(p,1),kk=wrap(k-olb::descriptors::c<D>(p,2),nz);if(py)jj=wrap(jj,ny);
        if(g.get({ii,jj,kk})==2)value+=cell[olb::descriptors::opposite<D>(p)];
        else {auto upstream=b.get(ii,jj,kk);value+=upstream[p];}
      }predicted[index(i,j,k)]=value;
    }
  }
  template<typename BLOCK,typename GEO>
  void afterStream(BLOCK& b,GEO& g,int step){
    T balance=0,reflection=0;std::size_t count=0;
    for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k)if(g.get({i,j,k})==1){
      auto cell=b.get(i,j,k);balance=std::max(balance,std::abs(olb::lbm<D>::computeRho(cell)-predicted[index(i,j,k)]));
      for(int p=1;p<D::q;++p){int ii=wrap(i+olb::descriptors::c<D>(p,0),nx),jj=j+olb::descriptors::c<D>(p,1),kk=wrap(k+olb::descriptors::c<D>(p,2),nz);if(py)jj=wrap(jj,ny);
        if(g.get({ii,jj,kk})==2){auto solid=b.get(ii,jj,kk);reflection=std::max(reflection,std::abs(cell[olb::descriptors::opposite<D>(p)]-solid[p]));++count;}
      }
    }
    std::ofstream out("halfway_audit_"+std::to_string(step)+".json");out<<std::setprecision(17)<<"{\"step\":"<<step<<",\"reflected_links\":"<<count<<",\"max_local_fluid_balance_residual\":"<<balance<<",\"max_reflection_pair_residual\":"<<reflection<<"}\n";
    if(balance>1e-11||reflection>1e-11)throw std::runtime_error("Stock half-way population diagnostic failed.");
  }
};
}
#endif
