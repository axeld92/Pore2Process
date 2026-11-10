using Test
include("AdsorptionColumn.jl")
using .AdsorptionColumn

@testset "Adsorption SI, limits and exact conservation" begin
    @test equilibrium_loading(.05,.04,100.)≈.04*5/6
    @test ldf_rate(1e-12,.002)≈1.5e-5
    @test pressure_load(1e-10,.001,.1,.001,.01,.75).pressure_drop_Pa≈1000.
    base=(;length_m=.1,area_m2=.001,porosity=.5,particle_density_kg_m3=1000.,
        superficial_velocity_m_s=.001,inlet_concentration_kg_m3=.05,t_end_s=200.,cells=80,dt_s=.1)
    v=solve_column(;base...,capacity_kg_kg=0.)
    @test v.breakthrough_time_s<50<v.sharp_front_capacity_time_s+1e-10
    @test v.outlet_concentration_kg_m3[end]≈.05 atol=1e-10
    @test v.normalized_mass_residual<1e-10
    @test all(>=(0),v.final_concentration_kg_m3)
    @test all(iszero,v.final_loading_kg_kg)
    @test v.inventory_kg≈v.adsorbed_inventory_kg+v.liquid_inventory_kg
    @test v.solid_utilization_at_breakthrough===nothing
    frozen=solve_column(;base...,ldf_s_inv=0.)
    @test frozen.outlet_concentration_kg_m3≈v.outlet_concentration_kg_m3
    slow=solve_column(;base...,ldf_s_inv=1e-4)
    @test slow.outlet_concentration_kg_m3[end]<v.outlet_concentration_kg_m3[end]
    @test slow.normalized_mass_residual<1e-10
    # Closed, well-mixed batch with linear isotherm: independent exact decay.
    K=.001;k=.01;eps=.5;rb=500.;cin=.05;t=50.
    analytical=cin*(eps/(eps+rb*K)+rb*K/(eps+rb*K)*exp(-k*(1+rb*K/eps)*t))
    errors=Float64[]
    for dt in (.5,.25,.125)
        a=solve_column(;base...,superficial_velocity_m_s=0.,t_end_s=t,dt_s=dt,
            initial_concentration_kg_m3=cin,linear_partition_m3_kg=K,ldf_s_inv=k)
        push!(errors,abs(a.final_concentration_kg_m3[1]-analytical))
        @test a.normalized_mass_residual<1e-10
        @test a.breakthrough_time_s===nothing
        @test maximum(a.final_concentration_kg_m3)-minimum(a.final_concentration_kg_m3)<1e-12
    end
    @test 1.9<errors[1]/errors[2]<2.1
    @test 1.9<errors[2]/errors[3]<2.1
    # Nonlinear closed-batch equilibrium follows an independent quadratic mass balance.
    m0=.5*.05;B=.5+500*.04169*94-m0*94;coef=.5*94
    root=2m0/(B+sqrt(B^2+4coef*m0))
    nonlinear=solve_column(;base...,superficial_velocity_m_s=0.,initial_concentration_kg_m3=.05,
        t_end_s=2000.,dt_s=1.,ldf_s_inv=.01)
    @test nonlinear.final_concentration_kg_m3[1]≈root rtol=1e-6
    @test nonlinear.final_loading_kg_kg[1]≈equilibrium_loading(root,.04169,94.) rtol=1e-6
    @test nonlinear.normalized_mass_residual<1e-9
    stiff=solve_column(superficial_velocity_m_s=.00026,ldf_s_inv=.0012,cells=640,dt_s=250.,t_end_s=5000.)
    @test all(>=(0),stiff.final_concentration_kg_m3)
    @test all(>=(0),stiff.final_loading_kg_kg)
    @test stiff.normalized_mass_residual<1e-9
    @test stiff.breakthrough_time_s===nothing
    # Linear local-equilibrium capacity gives an exact retarded plug-front time.
    for n in (40,80,160)
        a=solve_column(;base...,linear_partition_m3_kg=.001,equilibrium=true,cells=n,t_end_s=160.,dt_s=.1,breakthrough_fraction=.5)
        @test abs(a.breakthrough_time_s/100-1)<.03
        @test a.normalized_mass_residual<1e-10
    end
    for kwargs in ((;porosity=1.),(;cells=1),(;dt_s=NaN),(;affinity_m3_kg=-1.),(;inlet_concentration_kg_m3=Inf),
            (;ldf_s_inv=-1.),(;particle_density_kg_m3=0.),(;breakthrough_fraction=1.))
        @test_throws ArgumentError solve_column(;kwargs...)
    end
    @test_throws ArgumentError pressure_load(0.,.001,1.,.1,.1,.8)
end
