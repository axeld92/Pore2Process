module Pore2ProcessDemo

using LinearAlgebra, SparseArrays
export run_case, predict_column_flow, infer_permeability, capillary_reference, predict_column_pressure_drop
export filtration_load, permeance_LMH_bar_to_SI, resistance_from_permeance
export support_permeability, fit_support_multiplier, permeability_errors
export series_impedance_fit, series_flux, impedance_load
export graded_pore_reference, graded_support_load
export once_through_filtration
export crossflow_recycle_load
export crossflow_benchmark_point, critical_flux_bracket, crossflow_solute_balance
export duct_stokes_reference, star_membrane_load, star_feed_balance

function positive(name, x)
    x isa Real && !(x isa Bool) && isfinite(x) && x > 0 ||
        throw(ArgumentError("$name must be finite and greater than zero (SI units)."))
    y = Float64(x)
    isfinite(y) && y > 0 || throw(ArgumentError("$name is outside the supported Float64 range."))
    y
end

function finite(name, x)
    x isa Real && !(x isa Bool) && isfinite(x) || throw(ArgumentError("$name must be a finite number."))
    y = Float64(x)
    isfinite(y) || throw(ArgumentError("$name is outside the supported Float64 range."))
    y
end

function dimension(name, x, minimum)
    x isa Integer && !(x isa Bool) && minimum <= x <= typemax(Int) ||
        throw(ArgumentError("$name must be an integer >= $minimum."))
    Int(x)
end

"""
    run_case(nx=6, ny=4, nz=4, spacing_m=1e-3, radius_m=50e-6,
             viscosity_Pa_s=1e-3, p_in_Pa=100.0, p_out_Pa=0.0)

Solve the sparse capillary-network equations in SI units. Logical (i,j,k) has
one-based node ID 1+i+nx*(j+ny*k): x varies fastest. At each node, emit existing
positive-x, then positive-y, then positive-z edges. A positive edge flow travels
from its first to its second endpoint. Each call owns all its returned arrays.
"""
function run_case(nx=6, ny=4, nz=4, spacing_m=1e-3, radius_m=50e-6,
                  viscosity_Pa_s=1e-3, p_in_Pa=100.0, p_out_Pa=0.0)
    started = time_ns()
    nx, ny, nz = dimension("nx", nx, 2), dimension("ny", ny, 1), dimension("nz", nz, 1)
    a, r, mu = positive("spacing_m", spacing_m), positive("radius_m", radius_m), positive("viscosity_Pa_s", viscosity_Pa_s)
    pin, pout = finite("p_in_Pa", p_in_Pa), finite("p_out_Pa", p_out_Pa)
    2r < a || throw(ArgumentError("Capillary diameter must be smaller than spacing: require 2*radius_m < spacing_m."))
    dp = finite("pressure difference", pin-pout)
    N = Base.checked_mul(Base.checked_mul(nx, ny), nz)
    M = Base.checked_add(Base.checked_add(Base.checked_mul(nx-1, ny*nz), Base.checked_mul(nx, (ny-1)*nz)), Base.checked_mul(nx*ny, nz-1))
    L, A = positive("sample length", (nx-1)*a), positive("sample area", ny*nz*a^2)
    g = positive("throat conductance (adjust radius/viscosity/spacing to avoid overflow or underflow)", pi*r^4/(8mu*a))
    coordinates = Matrix{Float64}(undef, N, 3)
    edges = Matrix{Int}(undef, M, 2)
    node(i,j,k) = 1+i+nx*(j+ny*k)
    e = 0
    for k in 0:nz-1, j in 0:ny-1, i in 0:nx-1
        u = node(i,j,k)
        coordinates[u,:] = [i*a, (j+0.5)*a, (k+0.5)*a]
        for (exists, v) in ((i < nx-1, u+1), (j < ny-1, u+nx), (k < nz-1, u+nx*ny))
            if exists
                e += 1
                edges[e,:] = [u,v]
            end
        end
    end
    inlet = collect(1:nx:N)
    outlet = collect(nx:nx:N)
    boundary = vcat(inlet, outlet)
    interior = [u for u in 1:N if mod(u-1,nx) != 0 && mod(u-1,nx) != nx-1]
    # Assemble K from edge contributions; no analytical solution is used here.
    rows, cols, vals = Int[], Int[], Float64[]
    sizehint!(rows,4M); sizehint!(cols,4M); sizehint!(vals,4M)
    for t in 1:M
        u,v = edges[t,1], edges[t,2]
        append!(rows, (u,v,u,v)); append!(cols, (u,v,v,u)); append!(vals, (g,g,-g,-g))
    end
    K = sparse(rows, cols, vals, N, N)
    pressures = fill(pout,N)
    pressures[inlet] .= pin
    solve_started = time_ns()
    # nx=2 has no interior unknowns. At dp=0 set the exact constant state.
    if !isempty(interior) && dp != 0
        pressures[interior] = K[interior,interior] \ (-K[interior,boundary]*pressures[boundary])
    end
    solve_s = (time_ns()-solve_started)/1e9
    flows = [g*(pressures[edges[t,1]]-pressures[edges[t,2]]) for t in 1:M]
    net_outflows = zeros(N)
    for t in 1:M
        net_outflows[edges[t,1]] += flows[t]
        net_outflows[edges[t,2]] -= flows[t]
    end
    Qin, Qout = sum(net_outflows[inlet]), -sum(net_outflows[outlet])
    imbalance = isempty(interior) ? 0.0 : maximum(abs, net_outflows[interior])
    mismatch = abs(Qin-Qout)
    k_eff = dp == 0 ? nothing : mu*L*Qin/(A*dp)
    all(isfinite, pressures) && all(isfinite, flows) && all(isfinite, net_outflows) &&
        isfinite(Qin) && isfinite(Qout) && (isnothing(k_eff) || isfinite(k_eff)) ||
        throw(ArgumentError("Non-finite numerical result; use less extreme inputs."))
    scale = max(abs(Qin),abs(Qout))
    (; coordinates_m=coordinates, edges, pressures_Pa=pressures, flows_m3_s=flows,
       radii_m=fill(r,M), lengths_m=fill(a,M), conductances_m3_Pa_s=fill(g,M),
       net_outflows_m3_s=net_outflows, inlet_nodes=inlet, outlet_nodes=outlet, interior_nodes=interior,
       Q_in_m3_s=Qin, Q_out_m3_s=Qout, length_m=L, area_m2=A,
       permeability_m2=k_eff, max_interior_imbalance_m3_s=imbalance, flow_mismatch_m3_s=mismatch,
       normalized_conservation_error=scale == 0 ? 0.0 : max(imbalance,mismatch)/scale,
       linear_solve_s=solve_s, total_s=(time_ns()-started)/1e9)
end

"""Darcy column demonstration. `nothing` permeability propagates as unavailable."""
function predict_column_flow(permeability_m2, viscosity_Pa_s, length_m, area_m2, delta_p_Pa)
    mu, L, A = positive("column viscosity_Pa_s",viscosity_Pa_s), positive("column length_m",length_m), positive("column area_m2",area_m2)
    dp = finite("column delta_p_Pa",delta_p_Pa)
    isnothing(permeability_m2) && return nothing
    k = finite("permeability_m2",permeability_m2)
    k >= 0 || throw(ArgumentError("permeability_m2 must be nonnegative."))
    finite("column flow", k*A*dp/(mu*L))
end

include("Hydraulics.jl")
include("Filtration.jl")
include("Calibration.jl")
include("SeriesCalibration.jl")
include("GradedSupport.jl")
include("ProcessFiltration.jl")
include("CrossflowRecycle.jl")
include("CrossflowBenchmark.jl")
include("StarChannel.jl")

end
