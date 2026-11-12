// SPDX-License-Identifier: GPL-2.0-or-later
// Initial-data function only. Geometry/projection uses stock OpenLB facilities.
#ifndef P2P_COMPATIBLE_PROFILE_H
#define P2P_COMPATIBLE_PROFILE_H
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

class CompatibleInitialProfile {
  using T=double;
  T _side,_radius,_mean,_amplitude,_delta;
  olb::PolynomialStartScale<T,T> _blend{1.,1.};
  std::vector<olb::IndicatorSphere3D<T>> _spheres;
  T original(T x)const{return _mean+_amplitude*std::sin(2*M_PI*x/_side);}
public:
  CompatibleInitialProfile(T side,T mean,T amplitude)
    :_side(side),_radius(side/std::cbrt(2*M_PI/(3*.4))/2),
     _mean(mean),_amplitude(amplitude),_delta(.4*(side/std::sqrt(2.)-2*_radius)) {
    if(_delta<=0)throw std::invalid_argument("Compatibility buffer needs separated FCC spheres.");
    for(int i=-1;i<=3;++i)for(int j=-1;j<=3;++j)for(int k=-1;k<=3;++k){
      if((i+j+k)%2)continue;
      _spheres.emplace_back(olb::Vector<T,3>{i*side/2,j*side/2,k*side/2},_radius);
    }
  }
  T buffer()const{return _delta;}
  T operator()(olb::Vector<T,3> x){
    for(int dim=0;dim<3;++dim){x[dim]=std::fmod(x[dim],_side);if(x[dim]<0)x[dim]+=_side;}
    T distance=std::numeric_limits<T>::infinity();int nearest=-1;
    for(std::size_t i=0;i<_spheres.size();++i){
      T d=_spheres[i].signedDistance(x);
      if(d<distance){distance=d;nearest=i;}
    }
    if(distance>=_delta)return original(x[0]);
    auto& sphere=_spheres[nearest];auto radial=x-sphere.getCenter();
    if(olb::norm(radial)==0)return original(x[0]); // inactive deep-solid centre only
    auto surface=sphere.getCenter()+olb::normalize(radial,sphere.getRadius());
    T t=std::clamp(distance/_delta,T{0},T{1});
    T out[1];_blend(out,&t);T w=out[0];
    return (1-w)*original(surface[0])+w*original(x[0]);
  }
  void writeChecks(){
    std::ofstream meta("initial_profile_metadata.json");
    meta<<std::setprecision(17)<<"{\"name\":\"surface-trace preserving normal extension with quintic blend\",\"mean_kg_m3\":"<<_mean
      <<",\"amplitude_kg_m3\":"<<_amplitude<<",\"buffer_m\":"<<_delta
      <<",\"minimum_sphere_gap_m\":"<<_side/std::sqrt(2.)-2*_radius
      <<",\"sphere_radius_m\":"<<_radius<<",\"side_m\":"<<_side
      <<",\"probe_normal_step_m\":"<<_delta*1e-4
      <<",\"geometry_facilities\":\"stock IndicatorSphere3D signedDistance/getCenter/getRadius and olb::normalize\"}\n";
    std::ofstream out("initial_profile_probes.csv");out<<std::setprecision(17);
    out<<"group_id,sphere_id,nx,ny,nz,normal_distance_m,x_m,y_m,z_m,C_kg_m3,periodic_x_C,periodic_y_C,periodic_z_C\n";
    std::vector<olb::Vector<T,3>> centres={{0,0,0},{0,_side/2,_side/2},{_side/2,0,_side/2},{_side/2,_side/2,0}};
    std::vector<olb::Vector<T,3>> directions={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
      {1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{-1,-1,1},{-1,1,-1},{1,-1,-1},{-1,-1,-1},
      {1,2,3},{-1,2,3},{1,-2,3},{1,2,-3}};
    int group=0;
    for(std::size_t b=0;b<centres.size();++b)for(auto v:directions){
      v=olb::normalize(v);T h=_delta*1e-4;
      for(T d:std::vector<T>{-h,0,h,_delta-h,_delta,_delta+h}){
        auto x=centres[b]+(_radius+d)*v;
        auto xx=x,xy=x,xz=x;xx[0]+=_side;xy[1]+=_side;xz[2]+=_side;
        out<<group<<','<<b<<','<<v[0]<<','<<v[1]<<','<<v[2]<<','<<d<<','<<x[0]<<','<<x[1]<<','<<x[2]<<','
          <<(*this)(x)<<','<<(*this)(xx)<<','<<(*this)(xy)<<','<<(*this)(xz)<<'\n';
      }
      ++group;
    }
  }
};
#endif
