"""Independent checks of the new prescribed-flow verification oracles."""
import numpy as np

from examples.adsorption.transport_verification import reference as r


def test_uniform_mode_reduces_to_analytical_translation_and_decay():
    times=(.041,.171)
    f=r.modes(32,12,times,ax=0.,b=0.)
    for a,t in zip(f,times):
        np.testing.assert_allclose(a,np.exp((-r.D*r.K**2-1j*r.K*r.U)*t),atol=1e-13,rtol=0)


def test_stretching_at_zero_diffusion_matches_backward_characteristics():
    # Independently integrate the full characteristic vector with RK4.
    y=np.linspace(.023,.969,13)*r.L;z=np.linspace(.019,.978,17)*r.L
    Y,Z=np.meshgrid(y,z,indexing='ij');end=np.stack((np.zeros_like(Y),Y,Z),axis=-1)
    start=end.copy();time=.046;dt=-time/600
    def flow(a):return r.velocity('stretch',a[...,0],a[...,1],a[...,2])
    for _ in range(600):
        a=flow(start);b=flow(start+.5*dt*a);c=flow(start+.5*dt*b);d=flow(start+dt*c)
        start+=dt*(a+2*b+2*c+d)/6
    expected=np.exp(1j*r.K*start[...,0])
    f=r.modes(128,256,(time,),diffusivity=0.)[0]
    actual=r.sample(f,y,z)
    assert np.max(abs(actual-expected))<1e-5
    assert np.max(abs(np.diff(actual,axis=0)))>.1
    assert np.max(abs(np.diff(actual,axis=1)))>.1


def test_disk_root_and_zero_normal_flux():
    a=r.ROOT_J2_PRIME
    q=np.linspace(.01,a,1000)
    derivative=r.bessel(1,q)-2*r.bessel(2,q)/q
    assert min(derivative[:-1])>0
    assert abs(derivative[-1])<1e-13
    theta=np.linspace(0,2*np.pi,31)
    def c(radius):return r.disk_profile(r.L/2+radius*np.cos(theta),r.L/2+radius*np.sin(theta),.091)
    h=1e-8
    assert np.max(abs(c(r.R+h)-c(r.R-h)))*r.R/(2*h*r.A)<1e-8


def test_disk_reference_satisfies_advection_diffusion_equation():
    x=r.L/2+np.array([.06,-.2,.1])*1e-3;y=r.L/2+np.array([.22,.1,-.25])*1e-3;t=.087
    h=1e-8;dt=1e-6
    p=r.disk_profile
    c=p(x,y,t);xp=p(x+h,y,t);xm=p(x-h,y,t);yp=p(x,y+h,t);ym=p(x,y-h,t)
    ut=r.velocity('disk',x,y,np.zeros_like(x))
    residual=(p(x,y,t+dt)-p(x,y,t-dt))/(2*dt)+ut[:,0]*(xp-xm)/(2*h)+ut[:,1]*(yp-ym)/(2*h)-r.D*(xp+xm+yp+ym-4*c)/h**2
    assert np.max(abs(residual))<2e-6
