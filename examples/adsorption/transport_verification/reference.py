"""Continuum references for the two prescribed, divergence-free flows."""
import numpy as np

L = .0013891725543141453
U = .012937295182166382
D = 1e-9
DT = 1.9017455796408345e-5
STEPS = 11293
R = .0004
A = .25
AX = .5
B = .35
K = 2*np.pi/L
ROOT_J2_PRIME = 3.0542369282271403
TIMES = (STEPS//2*DT, STEPS*DT)


def velocity(kind, x, y, z, speed=U):
    if kind=='disk':
        return np.stack((-speed/R*(y-L/2),speed/R*(x-L/2),np.zeros_like(x)),axis=-1)
    if kind=='uniform':
        return np.stack((np.full_like(x,speed),np.zeros_like(x),np.zeros_like(x)),axis=-1)
    if kind!='stretch':raise ValueError('Unknown prescribed flow')
    return np.stack((speed*(1+AX*np.cos(K*y)*np.cos(K*z)),
                     B*speed*np.sin(K*z), B*speed*np.sin(K*y)),axis=-1)


def modes(points, steps_per_window, times=TIMES, *, diffusivity=D, ax=AX, b=B):
    """One exact x Fourier mode; split the two transverse shears in y,z."""
    q=(np.arange(points)+.5)*L/points
    ky=2*np.pi*np.fft.fftfreq(points,d=L/points)
    sin=np.sin(K*q);cos=np.cos(K*q)
    f=np.ones((points,points),dtype=complex);out=[];start=0.
    for end in times:
        dt=(end-start)/steps_per_window
        diff=np.exp(-.5*diffusivity*dt*(K*K+ky[:,None]**2+ky[None,:]**2))
        phase=np.exp(-.5j*dt*K*U*(1+ax*cos[:,None]*cos[None,:]))
        dy=np.exp(-.5j*dt*ky[:,None]*b*U*sin[None,:])
        dz=np.exp(-1j*dt*b*U*sin[:,None]*ky[None,:])
        for _ in range(steps_per_window):
            f=np.fft.ifft2(diff*np.fft.fft2(f))*phase
            f=np.fft.ifft(dy*np.fft.fft(f,axis=0),axis=0)
            f=np.fft.ifft(dz*np.fft.fft(f,axis=1),axis=1)
            f=np.fft.ifft(dy*np.fft.fft(f,axis=0),axis=0)*phase
            f=np.fft.ifft2(diff*np.fft.fft2(f))
        out.append(f.copy());start=end
    return np.asarray(out)


def sample(f, y, z):
    """Sample the continuum Fourier interpolant at cell centres or off-grid."""
    n=f.shape[0];freq=np.fft.fftfreq(n)*n
    def matrix(q):return np.exp(2j*np.pi*np.outer(np.asarray(q)/L-.5/n,freq))
    return matrix(y)@(np.fft.fft2(f)/n**2)@matrix(z).T


def qualify_stretch():
    """Separately refine space and time, with a finite and recorded budget."""
    cache={};comparisons=[];tol=1e-4
    probe=(np.arange(161)+.5)*L/161
    def get(n,s):
        if (n,s) not in cache:cache[n,s]=modes(n,s)
        return cache[n,s]
    def compare(n0,s0,n1,s1,axis):
        f0,f1=get(n0,s0),get(n1,s1)
        # Max over x is exactly the modulus of the complex envelope error.
        errors=[float(np.max(abs(sample(a,probe,probe)-sample(b,probe,probe)))) for a,b in zip(f0,f1)]
        comparisons.append(dict(axis=axis,coarse=[n0,s0],fine=[n1,s1],maximum_over_amplitude=errors))
        return max(errors)<=tol
    for n in (128,256,512):
        if not compare(n//2,256,n,256,'space'):continue
        for s in (512,1024):
            if not compare(n,s//2,n,s,'time'):continue
            if not compare(n//2,s,n,s,'space_at_selected_time_step'):break
            selected=get(n,s)
            info=dict(points=n,steps_per_window=s,times_s=list(TIMES),probe_points=161,
                      tolerance=tol,comparisons=comparisons,passed=True)
            return selected,info
    raise ValueError(f'Independent stretching reference missed its convergence screen: {comparisons}')


def gradients(f):
    """RMS derivatives of C/A over the entire cube, and an inverse-gradient scale."""
    k=2*np.pi*np.fft.fftfreq(len(f),d=L/len(f));F=np.fft.fft2(f)
    dy=np.fft.ifft2(1j*k[:,None]*F);dz=np.fft.ifft2(1j*k[None,:]*F)
    norms=[abs(K*f),abs(dy),abs(dz)]
    return dict(RMS_per_m=[float(np.sqrt(.5*np.mean(g*g))) for g in norms],
                inverse_gradient_lower_scale_m=float(1/np.max(np.sqrt(sum(g*g for g in norms)))))


def bessel(order, arg):
    """J0/J1/J2 power series on the small interval used by the disk oracle."""
    x=np.asarray(arg,dtype=float)
    if order not in (0,1,2) or np.any(abs(x)>8):raise ValueError('Outside disk oracle range')
    term=(x/2)**order/(2 if order==2 else 1);total=term.copy()
    for k in range(1,28):
        term=term*(-.25*x*x/(k*(k+order)));total+=term
    return total


def disk_profile(x, y, time, speed=U, amplitude=A):
    """Exact rigid rotation and diffusion of the first m=2 Neumann mode."""
    x=np.asarray(x)-L/2;y=np.asarray(y)-L/2
    arg=ROOT_J2_PRIME*np.hypot(x,y)/R
    theta=np.arctan2(y,x)-speed/R*time
    decay=np.exp(-D*(ROOT_J2_PRIME/R)**2*time)
    return .5+amplitude*decay*bessel(2,arg)/bessel(2,ROOT_J2_PRIME)*np.cos(2*theta)
