"""Independent high-Pe circular-tube first-passage moments.

This is a radial finite-volume representation of the analytic developed
Poiseuille velocity, not an OpenLB solver or a porous specimen. The axial
coordinate is the evolution variable. Flux-weighted injection makes the
radial distance process reversible with invariant flux weights. Molecular
diffusion is radial; axial molecular diffusion is neglected at high Pe.
"""
module TubeMoments

using LinearAlgebra

struct Spectrum
    radial_cells::Int
    weight::Vector{Float64}
    relative_velocity::Vector{Float64}
    mean_factor::Float64
    rates::Vector{Float64}
    covariance_weights::Vector{Float64}
    slope_dimensionless::Float64
end

function spectrum(n::Int; profile::Symbol=:parabolic)
    n >= 4 || throw(ArgumentError("At least four radial annuli are required"))
    edge=collect(range(0.0,1.0,length=n+1))
    area=diff(edge.^2)  # fraction of tube cross-sectional area
    # Exact area-average of 2*(1-r^2/R^2) over each annulus.
    velocity=profile==:parabolic ? 2 .- edge[1:end-1].^2 .- edge[2:end].^2 :
             profile==:plug ? ones(n) : throw(ArgumentError("Unknown velocity profile"))
    all(>(0),velocity) || error("Nonpositive averaged axial speed")
    w=area.*velocity
    @assert isapprox(sum(w),1;rtol=1e-13)
    # Dimensionless distance generator: w_i G_ij = conductance_ij = w_j G_ji.
    # Distance unit is ubar*R^2/D, not physical time.
    face=2n.*edge[2:end-1]
    diag=(vcat(0.,face)+vcat(face,0.))./w
    off=-face./sqrt.(w[1:end-1].*w[2:end])
    k=SymTridiagonal(diag,off)
    rootw=sqrt.(w)
    @assert norm(k*rootw,Inf)<=1e-11*maximum(diag)
    e=eigen(k)
    @assert abs(e.values[1])<=1e-11*maximum(e.values)
    rate=e.values[2:end]
    all(>(0),rate) || error("Nonpositive decaying mode")
    inverse_velocity=1 ./ velocity
    mean_factor=dot(w,inverse_velocity)
    coeff=e.vectors[:,2:end]'*(rootw.*(inverse_velocity .- mean_factor))
    weights=coeff.^2
    slope=2sum(weights./rate)
    Spectrum(n,w,velocity,mean_factor,rate,weights,slope)
end

# ∫_0^z (z-x)e^(-κx)dx = z² (κz + expm1(-κz))/(κz)².
# Small-x expansion prevents cancellation for short tubes/low modes.
@inline function integration_factor(x::Float64)
    x<1e-3 ? .5-x/6+x^2/24-x^3/120 : (x+expm1(-x))/x^2
end

function moments(spec::Spectrum,nominal_mean_s::Real,radius_m::Real,diffusivity_m2_s::Real)
    all(isfinite,(nominal_mean_s,radius_m,diffusivity_m2_s)) || throw(ArgumentError("Nonfinite tube input"))
    nominal_mean_s>=0 && radius_m>0 && diffusivity_m2_s>0 || throw(ArgumentError("Invalid tube inputs"))
    radial_scale_s=Float64(radius_m^2/diffusivity_m2_s)
    z=Float64(nominal_mean_s/radial_scale_s)
    # Time change: T = integral of inverse relative velocity over distance.
    # Var(T) = 2 * integral of the stationary covariance with triangular weight.
    variance=2*Float64(nominal_mean_s)^2*sum(spec.covariance_weights.*integration_factor.(spec.rates.*z))
    (mean_s=Float64(nominal_mean_s)*spec.mean_factor,variance_s2=variance,
     asymptotic_slope_s=spec.slope_dimensionless*radial_scale_s,
     continuum_slope_s=radial_scale_s/24)
end

function main(args)
    length(args)==4 || error("Use radial_cells_csv radius_m D_csv nominal_mean_times_csv")
    ns=parse.(Int,split(args[1],','));radius=parse(Float64,args[2])
    ds=parse.(Float64,split(args[3],','));times=parse.(Float64,split(args[4],','))
    BLAS.set_num_threads(1)
    println("radial_cells,diffusivity_m2_s,nominal_mean_s,mean_s,variance_s2,asymptotic_slope_s,continuum_slope_s")
    for n in ns
        e=spectrum(n)
        for d in ds,t in times
            m=moments(e,t,radius,d)
            println(join((n,d,t,m.mean_s,m.variance_s2,m.asymptotic_slope_s,m.continuum_slope_s),','))
        end
    end
end

end # module

if abspath(PROGRAM_FILE)==@__FILE__
    TubeMoments.main(ARGS)
end
