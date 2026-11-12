module RadialReference

using LinearAlgebra

const ALPHA = 3.0542369282271403

function j2(x)
    term = x^2 / 8
    total = term
    for k in 1:26
        term *= -x^2 / (4k * (k + 2))
        total += term
    end
    return total
end

initial(s) = j2(ALPHA * s) / j2(ALPHA)

function solve(n, substeps, times; diffusivity, radius, omega, shear=true)
    h = 1 / n
    edges = collect(range(0.0, 1.0; length=n+1))
    s = (edges[1:end-1] .+ edges[2:end]) ./ 2
    weights = diff(edges.^6) ./ 6
    kappa = diffusivity / radius^2
    # f=s²g removes the polar singularity. Radial diffusion of g has the
    # flux s⁵g′; the exact smooth wall condition is g′(1)=-2g(1).
    conductance = kappa .* edges[2:end-1].^5 ./ h
    lower = ComplexF64.(conductance ./ weights[2:end])
    upper = ComplexF64.(conductance ./ weights[1:end-1])
    omega_s = shear ? omega .* (1 .- s.^2) : fill(omega,n)
    diagonal = ComplexF64.(-2im .* omega_s)
    diagonal[1:end-1] .-= upper
    diagonal[2:end] .-= lower
    diagonal[end] -= 2kappa / ((1+h)*weights[end])
    operator = Tridiagonal(lower, diagonal, upper)
    g = ComplexF64.(initial.(s) ./ s.^2)
    scratch = similar(g)
    out = Matrix{ComplexF64}(undef,n,length(times))
    previous = 0.0
    for (i,t) in enumerate(times)
        dt = (t-previous)/substeps
        dt > 0 || error("Reference observation times must increase")
        left = lu(I - (dt/2)*operator)
        right = I + (dt/2)*operator
        for _ in 1:substeps
            mul!(scratch,right,g)
            ldiv!(left,scratch)
            g,scratch = scratch,g
        end
        out[:,i] = g .* s.^2
        previous = t
    end
    return s,out
end

function main(args)
    length(args)==8 || error("Use output.csv n steps t_mid t_final diffusivity radius omega")
    path,n,steps,t1,t2,d,r,omega=args
    s,f=solve(parse(Int,n),parse(Int,steps),(parse(Float64,t1),parse(Float64,t2));
              diffusivity=parse(Float64,d),radius=parse(Float64,r),omega=parse(Float64,omega))
    open(path,"w") do io
        println(io,"r_over_R,middle_real,middle_imag,final_real,final_imag")
        for i in eachindex(s)
            println(io,join((s[i],real(f[i,1]),imag(f[i,1]),real(f[i,2]),imag(f[i,2])),","))
        end
    end
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    RadialReference.main(ARGS)
end
