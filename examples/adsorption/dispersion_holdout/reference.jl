using LinearAlgebra
include(joinpath(@__DIR__,"../../../julia/adsorption/ColumnBoundaryReference.jl"))
using .ColumnBoundaryReference: danckwerts_moments

function fit_dispersion(length_m,velocity,variance,pulse_s)
    target = variance-pulse_s^2/12
    tau = length_m/velocity
    0 < target < tau^2 || error("Variance outside the scalar dispersion model range")
    variance_at(d) = danckwerts_moments(length_m=length_m,velocity_m_s=velocity,
                                       dispersion_m2_s=d).variance_s2
    lo,hi = 0.,velocity*length_m
    while variance_at(hi)<target
        hi *= 2
    end
    for _ in 1:64
        mid = (lo+hi)/2
        if variance_at(mid)<target
            lo = mid
        else
            hi = mid
        end
    end
    (lo+hi)/2
end

"""Conservative Crank–Nicolson pulse response with total-flux inlet."""
function predict(path,length_m,volume,flow,d,pulse,dt,steps,cells,substeps)
    area = volume/length_m
    v = flow/area
    h = length_m/cells
    delta = dt/substeps
    r = d/h^2
    a = v/(2h)
    v*h/d<=1 || error("Refine reference grid")
    diagonal = fill(-2r,cells)
    diagonal[1] = diagonal[end] = -r-a
    A = Tridiagonal(fill(r+a,cells-1),diagonal,fill(r-a,cells-1))
    F = lu(I-delta/2*A)
    G = I+delta/2*A
    c = zeros(cells)
    rhs = similar(c)
    removed = m1 = m2 = 0.
    dose = flow*pulse
    pulse_steps = round(Int,pulse/delta)
    isapprox(pulse_steps*delta,pulse;atol=1e-15,rtol=0) || error("Pulse event must match a time step")
    open(path,"w") do io
        println(io,"time_s,liquid_kg,min_kg_m3,max_kg_m3,injected_kg,recovered_kg,first_moment_kg_s,second_moment_kg_s2,interval_exit_kg")
        println(io,"0,0,0,0,0,0,0,0,0")
        for k in 1:steps
            before = removed
            for i in 1:substeps
                step = (k-1)*substeps+i
                t = (step-1)*delta
                f0 = flow*c[end]
                mul!(rhs,G,c)
                if step<=pulse_steps
                    rhs[1] += delta*v/h
                end
                ldiv!(F,rhs)
                c,rhs = rhs,c
                f1 = flow*c[end]
                removed += delta*(f0+f1)/2
                m1 += delta*((t/2+delta/6)*f0+(t/2+delta/3)*f1)
                m2 += delta*((t^2/2+t*delta/3+delta^2/12)*f0+
                             (t^2/2+2t*delta/3+delta^2/4)*f1)
            end
            time = k*dt
            entered = flow*min(time,pulse)
            liquid = area*h*sum(c)
            abs(entered-liquid-removed)<=1e-9*dose || error("Reference mass did not close")
            println(io,join((time,liquid,minimum(c),maximum(c),entered,removed,m1,m2,removed-before),','))
        end
    end
end

function driver(args)
    output = args[1]
    L,V,Q,variance,pulse,dt = parse.(Float64,args[2:7])
    short,long = parse.(Int,args[8:9])
    d = fit_dispersion(L,Q*L/V,variance,pulse)
    open(joinpath(output,"coefficient.txt"),"w") do io
        println(io,d)
    end
    for (n,steps) in ((1,short),(2,long))
        predict(joinpath(output,"model_$(n)_coarse.csv"),n*L,n*V,Q,d,pulse,dt,steps,200*n,2)
        predict(joinpath(output,"model_$(n).csv"),n*L,n*V,Q,d,pulse,dt,steps,400*n,4)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    driver(ARGS)
end
