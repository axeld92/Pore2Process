using Test
include("FilmParticle.jl")
using .FilmParticle

@testset "Separate finite-bath endpoint" begin
    c0=.2;V=.00025;m=.001;qm=.1115;b=139.;beta=m/V
    B=1+beta*qm*b-c0*b
    exact=2c0/(B+sqrt(B^2+4b*c0))
    e=bath_equilibrium()
    @test e.concentration_kg_m3≈exact rtol=1e-12
    @test e.loading_kg_kg≈(c0-exact)/beta rtol=1e-12
    @test abs(e.mass_residual_kg)<1e-18
    @test equilibrium(.15)≈qm*b*.15/(1+b*.15)
    h=initial_film_rate(particle_diameter_m=.0008,film_m_s=1e-5)
    small=initial_film_rate(particle_diameter_m=.0004,film_m_s=2e-5)
    @test small.initial_rate_kg_kg_s≈4h.initial_rate_kg_kg_s
end

@testset "Fixed-surface spherical diffusion and shell convergence" begin
    R=.0004;D=1e-12;t=20000.;K=.001;C=.2
    exact=K*C*(1-6/pi^2*sum(exp(-n^2*pi^2*D*t/R^2)/n^2 for n in 1:100))
    errors=Float64[]
    for nr in (8,16,32)
        r=solve_fixed_bulk(concentration_kg_m3=C,particle_diameter_m=2R,
            diffusivity_m2_s=D,film_m_s=1.,linear_partition_m3_kg=K,
            radial_cells=nr,t_end_s=t,dt_s=5.)
        push!(errors,abs(last(r.mean_loading_kg_kg)-exact))
        @test r.normalized_mass_residual<1e-10
        @test all(>=(0),r.final_shell_loading_kg_kg)
        @test sum(r.radial_volume_weights)≈1
        @test size(r.final_shell_loading_kg_kg)==(1,nr)
    end
    @test errors[end]/(K*C)<.001
    @test errors[3]<errors[1]
end

@testset "Finite bath, film limit and conservative sampling" begin
    p=(;particle_diameter_m=.0008,particle_density_kg_m3=1000.,
        diffusivity_m2_s=1e-8,film_m_s=1e-8,linear_partition_m3_kg=.001,
        initial_concentration_kg_m3=.2,fluid_volume_m3=.00025,adsorbent_mass_kg=.001,
        radial_cells=12,t_end_s=7200.,dt_s=5.)
    r=solve_bath(;p...)
    h=3p.film_m_s/(p.particle_density_kg_m3*(p.particle_diameter_m/2))
    beta=p.adsorbent_mass_kg/p.fluid_volume_m3;K=p.linear_partition_m3_kg
    exact=.2/(beta+1/K)*(1-exp(-h*(beta+1/K)*p.t_end_s))
    @test last(r.mean_loading_kg_kg)≈exact rtol=.001
    @test last(r.outlet_concentration_kg_m3)≈.2-beta*last(r.mean_loading_kg_kg) atol=1e-10
    @test r.normalized_mass_residual<1e-8
    sampled=solve_bath(;p...,sampling_times_s=[3600.,7200.],sample_volume_m3=1e-6)
    @test last(sampled.fluid_volume_m3)≈.000248
    @test last(sampled.cumulative_sampled_solute_kg)>0
    @test sampled.normalized_mass_residual<1e-8
    nonlinear=solve_bath(t_end_s=2e6,dt_s=500.,diffusivity_m2_s=1e-10,film_m_s=1e-5)
    @test last(nonlinear.mean_loading_kg_kg)≈bath_equilibrium().loading_kg_kg rtol=1e-5
    @test nonlinear.normalized_mass_residual<1e-8
end

@testset "Column limits, conservation and independently owned runs" begin
    p=(;length_m=.1,area_m2=.001,porosity=.5,particle_density_kg_m3=1000.,
        superficial_velocity_m_s=.001,inlet_concentration_kg_m3=.05,
        t_end_s=200.,dt_s=.1,axial_cells=40,radial_cells=8)
    passive=solve_column(;p...,diffusivity_m2_s=0.)
    @test passive.normalized_mass_residual<1e-8
    @test last(passive.outlet_concentration_kg_m3)≈.05 atol=1e-10
    @test all(iszero,passive.final_shell_loading_kg_kg)
    blocked=solve_column(;p...,film_m_s=0.)
    @test blocked.outlet_concentration_kg_m3≈passive.outlet_concentration_kg_m3
    ads=solve_column(;p...,diffusivity_m2_s=1e-11,film_m_s=1e-5)
    @test last(ads.outlet_concentration_kg_m3)<last(passive.outlet_concentration_kg_m3)
    @test ads.normalized_mass_residual<1e-8
    @test ads.normalized_interface_flux_mismatch<1e-8
    @test ads.inventory_kg≈ads.liquid_inventory_kg+ads.adsorbed_inventory_kg
    @test all(>=(0),ads.final_concentration_kg_m3)
    @test size(ads.final_shell_loading_kg_kg)==(40,8)
    @test ads.physical_adsorbent_mass_kg≈.05
    @test ads.active_adsorbent_mass_kg≈.05
    half=solve_column(;p...,active_inventory_fraction=.5)
    @test half.physical_adsorbent_mass_kg==ads.physical_adsorbent_mass_kg
    @test half.active_adsorbent_mass_kg≈.025
    passive.final_shell_loading_kg_kg[1]=1.
    @test iszero(blocked.final_shell_loading_kg_kg[1])
    for invalid in ((;porosity=1.),(;particle_diameter_m=0.),(;diffusivity_m2_s=-1.),
            (;film_m_s=NaN),(;radial_cells=0),(;axial_cells=true),(;active_inventory_fraction=1.1))
        @test_throws ArgumentError solve_column(;invalid...)
    end
    @test_throws ArgumentError solve_bath(sampling_times_s=[2.,1.])
    @test_throws ArgumentError solve_bath(sampling_times_s=[1.],sample_volume_m3=.001)
end
