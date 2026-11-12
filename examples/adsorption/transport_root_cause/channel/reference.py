"""Independent continuum reference for one streamwise Fourier mode.

This is a verification oracle, not a second pore-scale production solver.
The divergence-free velocity is u=(U[1+a*cos(2*pi*y/H)],0,0). For
C=.5+A*Im[exp(i*k*x)*f(y,t)], the exact PDE reduces to
f_t=D*(f_yy-k*k*f)-i*k*u(y)*f, f(y,0)=1, periodic in y.
FFT diffusion and exact pointwise phase evolution use Strang splitting.
"""
import numpy as np


def mode(*, side_m, speed_m_s, shear_fraction, diffusivity_m2_s, duration_s,
         points=256, steps=2048):
    if points < 16 or steps < 1 or side_m <= 0 or duration_s < 0 or diffusivity_m2_s <= 0:
        raise ValueError('Positive physical scales, >=16 Fourier points and >=1 step required.')
    y=(np.arange(points)+.5)*side_m/points
    k=2*np.pi/side_m
    dt=duration_s/steps
    velocity=speed_m_s*(1+shear_fraction*np.cos(k*y))
    wave=2*np.pi*np.fft.fftfreq(points,d=side_m/points)
    diffusion_half=np.exp(-.5*diffusivity_m2_s*(wave*wave+k*k)*dt)
    phase=np.exp(-1j*k*velocity*dt)
    f=np.ones(points,dtype=complex)
    for _ in range(steps):
        f=np.fft.ifft(np.fft.fft(f)*diffusion_half)
        f*=phase
        f=np.fft.ifft(np.fft.fft(f)*diffusion_half)
    return f


def sample_mode(f, y_m, side_m):
    """Trigonometric interpolation of a periodic cell-centred Fourier grid."""
    count=len(f)
    frequencies=np.fft.fftfreq(count)*count
    coefficients=np.fft.fft(f)/count
    coordinates=np.asarray(y_m)/side_m-.5/count
    return np.exp(2j*np.pi*np.outer(coordinates,frequencies))@coefficients


def concentration(x_m,y_m, *, amplitude_kg_m3=.25, points=256, steps=2048, **physical):
    f=mode(points=points,steps=steps,**physical)
    values=sample_mode(f,y_m,physical['side_m'])
    return .5+amplitude_kg_m3*np.imag(np.exp(2j*np.pi*np.asarray(x_m)/physical['side_m'])*values)


def converged_reference(x_m,y_m, *, amplitude_kg_m3=.25, **physical):
    """Separate spatial and temporal checks; tolerance relative to amplitude."""
    coarse=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=128,steps=1024,**physical)
    spatial=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=256,steps=1024,**physical)
    fine=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=256,steps=2048,**physical)
    scale=max(amplitude_kg_m3,1e-30)
    errors=dict(spatial_max_over_amplitude=float(np.max(abs(coarse-spatial)))/scale,
                temporal_max_over_amplitude=float(np.max(abs(spatial-fine)))/scale,
                points=256,steps=2048,tolerance=1e-4)
    errors['passed']=max(errors['spatial_max_over_amplitude'],errors['temporal_max_over_amplitude'])<=errors['tolerance']
    if not errors['passed']:
        raise ValueError('Independent shear reference has not converged: '+str(errors))
    return fine,errors
