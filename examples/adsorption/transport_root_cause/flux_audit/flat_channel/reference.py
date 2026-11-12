"""Independent one-x-mode continuum oracle with Neumann channel walls.

Even extension of y in [0,H] to a period2H gives cosine/Neumann diffusion.
FFT diffusion and exact pointwise phase use Strang splitting. This is a
verification oracle, not a second production porous transport solver.
"""
import numpy as np

from ...channel.reference import sample_mode


def mode(*, height_m, mean_speed_m_s, diffusivity_m2_s, duration_s,
         points=512, steps=2048):
    if height_m<=0 or mean_speed_m_s<0 or diffusivity_m2_s<=0 or duration_s<0 or points<16 or steps<1:
        raise ValueError('Positive H/D, nonnegative U/time, >=16 points and >=1 step required')
    period=2*height_m
    y=(np.arange(2*points)+.5)*height_m/points
    reflected=np.minimum(y,period-y)/height_m
    velocity=6*mean_speed_m_s*reflected*(1-reflected)
    k=2*np.pi/height_m;dt=duration_s/steps
    wave=2*np.pi*np.fft.fftfreq(2*points,d=height_m/points)
    diffusion=np.exp(-.5*diffusivity_m2_s*(wave*wave+k*k)*dt)
    phase=np.exp(-1j*k*velocity*dt)
    f=np.ones(2*points,dtype=complex)
    for _ in range(steps):
        f=np.fft.ifft(np.fft.fft(f)*diffusion);f*=phase
        f=np.fft.ifft(np.fft.fft(f)*diffusion)
    return f


def concentration(x_m,y_m,*,amplitude_kg_m3=.25,points=512,steps=2048,**physical):
    f=mode(points=points,steps=steps,**physical)
    unique_y,inverse=np.unique(np.asarray(y_m),return_inverse=True)
    sampled=sample_mode(f,unique_y,2*physical['height_m'])[inverse]
    return .5+amplitude_kg_m3*np.imag(np.exp(2j*np.pi*np.asarray(x_m)/physical['height_m'])*sampled)


def boundary_checks(f,height_m):
    """Evaluate the exact derivative of the reference trigonometric polynomial."""
    points=len(f);frequencies=np.fft.fftfreq(points)*points
    coefficient=np.fft.fft(f)/points
    y=np.array([0.,height_m]);coordinates=y/(2*height_m)-.5/points
    derivative=np.exp(2j*np.pi*np.outer(coordinates,frequencies))@(coefficient*(2j*np.pi*frequencies/(2*height_m)))
    return dict(even_reflection_error=float(np.max(abs(f-f[::-1]))),
                max_wall_normal_mode_derivative_times_height=float(np.max(abs(derivative)))*height_m)


def converged_reference(x_m,y_m,*,amplitude_kg_m3=.25,**physical):
    coarse=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=1024,steps=2048,**physical)
    spatial=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=2048,steps=2048,**physical)
    fine=concentration(x_m,y_m,amplitude_kg_m3=amplitude_kg_m3,points=2048,steps=4096,**physical)
    scale=max(amplitude_kg_m3,1e-30)
    checks=boundary_checks(mode(points=2048,steps=4096,**physical),physical['height_m'])
    result=dict(spatial_max_over_amplitude=float(np.max(abs(coarse-spatial)))/scale,
                temporal_max_over_amplitude=float(np.max(abs(spatial-fine)))/scale,
                points_on_physical_height=2048,points_on_even_period=4096,steps=4096,tolerance=1e-4,
                **checks)
    result['passed']=max(result['spatial_max_over_amplitude'],result['temporal_max_over_amplitude'])<=1e-4 and checks['even_reflection_error']<=1e-11 and checks['max_wall_normal_mode_derivative_times_height']<=1e-9
    if not result['passed']:raise ValueError('Neumann continuum reference not qualified: '+str(result))
    return fine,result
