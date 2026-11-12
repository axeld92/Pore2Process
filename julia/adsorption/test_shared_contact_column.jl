using Test
include("FilmParticle.jl")
include("ExchangeColumn.jl")
include("SharedContactColumn.jl")

@testset "Shared carbon with separate liquid films" begin
    p = (;length_m=.10,area_m2=pi*.011^2/4,porosity=.427,superficial_velocity_m_s=.00768,
          inlet_concentration_kg_m3=.15,dt_s=.05,t_end_s=20.,radial_cells=12,
          diffusivity_m2_s=1.846e-12,film_m_s=4e-5)
    simple = FilmParticle.solve_column(;p...,axial_cells=4)
    one = SharedContactColumn.solve_column([1.],[1.],ones(1,1),[1.];p...,periods=2,subcells=2)
    @test one.outlet_concentration_kg_m3 ≈ simple.outlet_concentration_kg_m3 rtol=1e-8 atol=1e-13
    @test one.adsorbed_inventory_kg ≈ simple.adsorbed_inventory_kg rtol=1e-8
    mix = [.8 .2; .2 .8]
    w = [.5,.5]
    equal = SharedContactColumn.solve_column(w,[1.,1.],mix,w;p...,periods=2,subcells=2)
    @test equal.outlet_concentration_kg_m3 ≈ one.outlet_concentration_kg_m3 rtol=1e-11 atol=1e-13
    inert = SharedContactColumn.solve_column(w,[.5,1.5],mix,[0.,1.];p...,periods=2,subcells=2,diffusivity_m2_s=0.)
    old_inert = ExchangeColumn.solve_column(w,[.5,1.5],mix;p...,periods=2,subcells=2,diffusivity_m2_s=0.)
    @test inert.outlet_concentration_kg_m3 ≈ old_inert.outlet_concentration_kg_m3 atol=1e-14
    @test inert.cumulative_film_uptake_kg==0
    @test inert.cumulative_film_release_kg==0
    @test inert.adsorbed_inventory_kg[end]==0
    @test_throws ArgumentError SharedContactColumn.solve_column(w,[.5,1.5],mix,[.2,.9];p...,periods=2)

    # Independent coupled linear system for two paths and one radial shell.
    h,H,R,rho = .5,.4,.0004,740/(1-.427)
    f = [.3,.7]
    r = SharedContactColumn.solve_column(w,[.5,1.5],mix,f;p...,periods=1,subcells=1,
        dt_s=h,t_end_s=h,radial_cells=1,linear_partition_m3_kg=H)
    V = r.liquid_volume_m3.*w.*[.5,1.5]
    Q = p.area_m2*p.superficial_velocity_m_s.*w
    A = r.external_particle_area_m2
    K = p.film_m_s*A.*f
    ks = rho*p.diffusivity_m2_s/(R/2)
    b = h*A*ks/r.physical_adsorbent_mass_kg
    mat = [V[1]+h*Q[1]+h*K[1] 0 -h*K[1] 0 0;
           0 V[2]+h*Q[2]+h*K[2] 0 -h*K[2] 0;
           p.film_m_s 0 -(p.film_m_s+ks*H) 0 ks;
           0 p.film_m_s 0 -(p.film_m_s+ks*H) ks;
           0 0 -b*f[1]*H -b*f[2]*H 1+b]
    c1,c2,_,_,q = mat\[h*Q[1]*.15,h*Q[2]*.15,0,0,0]
    @test r.outlet_concentration_kg_m3[end] ≈ (c1+c2)/2 rtol=1e-10
    @test r.adsorbed_inventory_kg[end] ≈ q*r.physical_adsorbent_mass_kg rtol=1e-10
    @test r.normalized_mass_residual < 1e-10
    @test r.normalized_interface_flux_mismatch < 1e-10
    @test sum(r.contact_area_by_path_m2) ≈ A
end
