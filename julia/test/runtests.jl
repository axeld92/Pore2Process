using Test, Pore2ProcessDemo

include("hydraulics.jl")
include("filtration.jl")
include("calibration.jl")
include("series_calibration.jl")
include("graded_support.jl")
include("process_filtration.jl")
include("crossflow_recycle.jl")
include("crossflow_benchmark.jl")
include("star_channel.jl")

function reference(nx=6,ny=4,nz=4,a=1e-3,r=50e-6,mu=1e-3,pin=100.0,pout=0.0)
    s = run_case(nx,ny,nz,a,r,mu,pin,pout)
    L = (nx-1)*a
    Q = ny*nz*pi*r^4*(pin-pout)/(8mu*L)
    k = pi*r^4/(8a^2)
    @test size(s.coordinates_m) == (nx*ny*nz,3)
    @test size(s.edges) == ((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1),2)
    @test maximum(abs.(s.pressures_Pa .- (pin .- (pin-pout).*s.coordinates_m[:,1]./L))) <= 1e-9*max(abs(pin-pout),1)
    @test s.Q_in_m3_s ≈ Q rtol=1e-9
    @test s.Q_out_m3_s ≈ Q rtol=1e-9
    @test s.permeability_m2 ≈ k rtol=1e-9
    @test s.max_interior_imbalance_m3_s <= 1e-9*abs(Q)
    @test s.flow_mismatch_m3_s <= 1e-9*abs(Q)
    transverse = [t for t in axes(s.edges,1) if s.coordinates_m[s.edges[t,1],1] == s.coordinates_m[s.edges[t,2],1]]
    @test isempty(transverse) || maximum(abs,s.flows_m3_s[transverse]) <= 1e-9*abs(Q)
    @test all(1 .<= s.edges .<= nx*ny*nz)
    @test all(s.edges[:,1] .< s.edges[:,2])
    @test all(t -> count(!iszero, s.coordinates_m[s.edges[t,2],:] - s.coordinates_m[s.edges[t,1],:]) == 1, axes(s.edges,1))
    @test all(t -> sum(abs, s.coordinates_m[s.edges[t,2],:] - s.coordinates_m[s.edges[t,1],:]) ≈ a, axes(s.edges,1))
    @test length(unique(Tuple.(eachrow(s.edges)))) == size(s.edges,1)
    s
end

@testset "Scientific reference, geometry and scaling" begin
    s = reference()
    @test s.Q_in_m3_s ≈ 7.853981633974485e-10 rtol=1e-9
    @test s.conductances_m3_Pa_s[1] ≈ 2.4543692606170265e-12 rtol=1e-9
    @test s.length_m == 0.005
    @test s.area_m2 == 1.6e-5
    @test s.edges == run_case().edges
    @test s.coordinates_m == run_case().coordinates_m
    single = reference(2,1,1)
    @test single.max_interior_imbalance_m3_s == 0
    reference(2,4,3)
    reference(9,3,2)
    longer = reference(11)
    @test longer.Q_in_m3_s ≈ s.Q_in_m3_s/2 rtol=1e-9
    for (r,mu,pin,pout,qfactor,kfactor) in ((50e-6,1e-3,200.,0.,2.,1.), (50e-6,1e-3,117.,17.,1.,1.), (50e-6,1e-3,0.,100.,-1.,1.), (100e-6,1e-3,100.,0.,16.,16.), (50e-6,2e-3,100.,0.,0.5,1.))
        t=reference(6,4,4,1e-3,r,mu,pin,pout)
        @test t.Q_in_m3_s ≈ qfactor*s.Q_in_m3_s rtol=1e-9
        @test t.permeability_m2 ≈ kfactor*s.permeability_m2 rtol=1e-9
        if pin == 117.
            @test maximum(abs.(t.flows_m3_s-s.flows_m3_s)) < 1e-9*abs(s.Q_in_m3_s)
        end
    end
    @test predict_column_flow(s.permeability_m2,1e-3,0.02,1e-4,1000.) ≈ 1.227184630308513e-8 rtol=1e-9
    @test predict_column_flow(s.permeability_m2,1e-3,0.02,1e-4,-1000.) < 0
end

@testset "Zero drop and validation" begin
    for dims in ((6,4,4),(2,1,1))
        z=run_case(dims...,1e-3,50e-6,1e-3,73.,73.)
        @test all(==(73.),z.pressures_Pa)
        @test maximum(abs,z.flows_m3_s) < 1e-20
        @test z.Q_in_m3_s == z.Q_out_m3_s == 0
        @test isnothing(z.permeability_m2)
        @test isnothing(predict_column_flow(z.permeability_m2,1e-3,0.02,1e-4,1000.))
    end
    defaults=Any[6,4,4,1e-3,50e-6,1e-3,100.,0.]
    for (index,value) in ((1,1),(1,2.5),(1,true),(2,0),(3,-1),(4,0.),(5,-1.),(5,0.5e-3),(6,0.),(7,Inf),(8,NaN),(4,Inf),(5,NaN),(6,Inf),(1,NaN))
        bad=copy(defaults); bad[index]=value
        @test_throws ArgumentError run_case(bad...)
    end
    for args in ((-1.,1e-3,.02,1e-4,1000.),(1e-12,0.,.02,1e-4,1000.),(nothing,1e-3,0.,1e-4,1000.),(1e-12,1e-3,.02,NaN,1000.),(1e-12,1e-3,.02,1e-4,Inf))
        @test_throws ArgumentError predict_column_flow(args...)
    end
end
