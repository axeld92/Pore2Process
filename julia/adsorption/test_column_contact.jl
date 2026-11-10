using Test
include("FilmParticle.jl")
include("ColumnContact.jl")
using .ColumnContact

@testset "Independent column contact verification" begin
    # Derive the N-stage step from Q*(C[j-1]-C[j])=Vj*dC[j]/dt.
    base=(length_m=.01,area_m2=.001,porosity=.5,superficial_velocity_m_s=.001,
        inlet_concentration_kg_m3=.2,t_end_s=50.,dt_s=.002,radial_cells=4,film_m_s=0.)
    tau=base.porosity*base.length_m/base.superficial_velocity_m_s
    for n in (1,2,4,16)
        r=solve_contact(;base...,sections=n)
        exact=tracer_step.(r.time_s;sections=n,residence_time_s=tau)
        @test maximum(abs.(r.outlet_concentration_kg_m3/.2-exact))<.001
        @test r.normalized_mass_residual<1e-9
        @test r.section_liquid_volume_m3*n≈base.porosity*base.area_m2*base.length_m
        @test r.section_dry_mass_kg*n≈r.physical_adsorbent_mass_kg
        @test r.active_adsorbent_mass_kg==r.physical_adsorbent_mass_kg
        # Independent survival-integral identities: E[T]=int S dt,
        # E[T^2]=2*int t*S dt (no numerical differentiation of a step).
        t=collect(0.:.001:100.)
        S=1 .-tracer_step.(t;sections=n,residence_time_s=tau)
        trap(y)=.001*(sum(y)-(first(y)+last(y))/2)
        mean=trap(S);variance=2trap(t.*S)-mean^2
        @test isapprox(mean,tau;rtol=1e-6)
        @test isapprox(variance,tau^2/n;rtol=1e-5)
        if n==1
            @test maximum(abs.(exact-(1 .-exp.(-r.time_s/tau))))<1e-14
        end
    end
    # Fast-particle linear-isotherm limit is an analytically retarded RTD.
    # Total section storage is eps + rho_bulk*K, including liquid storage.
    p=(length_m=.01,area_m2=.001,porosity=.5,superficial_velocity_m_s=.001,
        inlet_concentration_kg_m3=.2,particle_density_kg_m3=1000.,
        linear_partition_m3_kg=.001,diffusivity_m2_s=.001,film_m_s=1000.,
        particle_diameter_m=.0008,t_end_s=80.,dt_s=.002,radial_cells=8,sections=4)
    r=solve_contact(;p...)
    retardation=1+(1-p.porosity)*p.particle_density_kg_m3*p.linear_partition_m3_kg/p.porosity
    exact=tracer_step.(r.time_s;sections=p.sections,residence_time_s=tau*retardation)
    @test maximum(abs.(r.outlet_concentration_kg_m3/.2-exact))<.001
    @test r.normalized_mass_residual<1e-8
    # Finite-rate nonlinear sorption: reuse the tested particle surface flux.
    r=solve_contact(sections=4,t_end_s=200.,dt_s=.2)
    @test r.normalized_mass_residual<1e-8
    @test r.normalized_interface_flux_mismatch<1e-8
    @test minimum(r.final_shell_loading_kg_kg)>=0
    @test maximum(r.final_shell_loading_kg_kg)<=.1115
    @test minimum(r.final_concentration_kg_m3)>=0
    @test r.physical_adsorbent_mass_kg≈740*(pi*.011^2/4)*.1
    r=solve_contact(sections=1,superficial_velocity_m_s=0.,t_end_s=1.,dt_s=.1)
    @test r.liquid_residence_time_s===nothing
    @test all(==(0.),r.inventory_kg)
    @test all(==(0.),r.cumulative_inlet_kg)
    a=solve_contact(;base...,sections=3)
    b=solve_contact(;merge(base,(film_m_s=1e-5,diffusivity_m2_s=0.))...,sections=3)
    @test a.outlet_concentration_kg_m3≈b.outlet_concentration_kg_m3
    for n in (0,129,1.5,true)
        @test_throws ArgumentError solve_contact(sections=n)
    end
    @test_throws ArgumentError solve_contact(porosity=1.)
    @test_throws ArgumentError solve_contact(superficial_velocity_m_s=NaN)
    @test_throws MethodError solve_contact(active_inventory_fraction=.5)
    @test tracer_step(0.;sections=64)==0.
    @test tracer_step(100.;sections=64)==1.
end
