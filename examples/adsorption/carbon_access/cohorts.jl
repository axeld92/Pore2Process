module AccessCohorts

using LinearAlgebra
using ..SharedContactColumn: particle_surface!

function setup(p,contact,ns,mass,area,stream,cin,fraction,ratio)
    0 <= fraction <= 1 && 0 <= ratio <= 1 ||
        throw(ArgumentError("Direct carbon fraction and slow access ratio must be between zero and one"))
    fractions = (Float64(fraction),1-Float64(fraction))
    film = (p.kf,p.kf*ratio)
    ks = p.rho*p.Ds/(p.R/(2p.nr))
    cap = p.qe(cin)
    active = ntuple(g->fractions[g]>0 && film[g]>0 && ks>0 && cap>0,2)
    ntu = ntuple(g->active[g] ? area*fractions[g]*film[g].*contact./stream :
                 zeros(size(contact)),2)
    beta = ntuple(g->active[g] ? (ks/film[g]).*(1 .+ ntu[g]./(1 .+ ntu[3-g])) :
                  zeros(size(contact)),2)
    n = size(contact,1)
    loading = ntuple(_->zeros(p.nr,ns),2)
    (;p,cin,cap,mass,area,contact,fractions,film,ks,active,ntu,beta,loading,
      interior=ntuple(_->similar(loading[1]),2),surface=(zeros(n),zeros(n)),
      qsurface=(zeros(n),zeros(n)),rhs=zeros(n),previous=zeros(n))
end

function diffuse!(bed,op)
    for g in 1:2
        bed.active[g] || continue
        copyto!(bed.interior[g],bed.loading[g])
        ldiv!(op.F,bed.interior[g])
    end
end

function surface!(bed,c,section,kind,g,op,qmax,b,linear)
    other = 3-g
    for i in eachindex(c)
        ntu = bed.ntu[other][i,kind]
        bed.rhs[i] = (clamp(c[i],0.,bed.cin)+ntu*bed.surface[other][i])/(1+ntu)
    end
    mean_surface,_ = particle_surface!(bed.surface[g],bed.qsurface[g],bed.rhs,
        view(bed.beta[g],:,kind),view(bed.contact,:,kind),bed.loading[g][end,section],
        bed.interior[g][end,section],op.b[end],bed.cap,qmax,b,linear,bed.p.qe,bed.p.slope)
    mean_surface
end

function solve_surfaces!(bed,c,section,kind,op,qmax,b,linear)
    for s in bed.surface
        fill!(s,0.)
    end
    if all(bed.active)
        for _ in 1:60
            copyto!(bed.previous,bed.surface[2])
            first = surface!(bed,c,section,kind,1,op,qmax,b,linear)
            second = surface!(bed,c,section,kind,2,op,qmax,b,linear)
            change = maximum(abs(bed.surface[2][i]-bed.previous[i]) for i in eachindex(c))
            change<=2e-12*max(bed.cin,1e-6) && return (first,second)
        end
        error("The two carbon surfaces did not converge to a common liquid concentration")
    end
    ntuple(g->bed.active[g] ? surface!(bed,c,section,kind,g,op,qmax,b,linear) : 0.,2)
end

"""Advance both intraparticle profiles with conservative fluxes from the same liquid."""
function contact!(bed,c,section,kind,h,op,qmax,b,linear)
    mean_surface = solve_surfaces!(bed,c,section,kind,op,qmax,b,linear)
    maxq = 0.
    adsorbed = ntuple(2) do g
        bed.active[g] || return 0.
        mean_loading = 0.
        for r in eachindex(op.w)
            q = bed.interior[g][r,section]+op.b[r]*mean_surface[g]
            q >= -1e-13 || error("Negative particle loading")
            bed.loading[g][r,section] = q
            maxq = max(maxq,q)
            mean_loading += op.w[r]*q
        end
        bed.mass*bed.fractions[g]*mean_loading
    end
    maxflux = uptake = release = 0.
    for i in eachindex(c)
        first,second = bed.ntu[1][i,kind],bed.ntu[2][i,kind]
        cout = (c[i]+first*bed.surface[1][i]+second*bed.surface[2][i])/(1+first+second)
        for g in 1:2
            bed.ntu[g][i,kind]>0 || continue
            film = bed.film[g]*(cout-bed.surface[g][i])
            radial = bed.ks*(bed.qsurface[g][i]-bed.loading[g][end,section])
            maxflux = max(maxflux,abs(film-radial))
            transfer = h*bed.area*bed.fractions[g]*bed.contact[i,kind]*film
            uptake += max(transfer,0.)
            release += max(-transfer,0.)
        end
        c[i] = cout
    end
    (;adsorbed,maxq,maxflux,uptake,release)
end

end
