module DiskFaceTransport

using LinearAlgebra
using SparseArrays

export sector_grid, transport

"""
Geometry and integrated face coefficients for an axially invariant disk.
The azimuthal flow is positive in the direction of increasing angle.
All radial and axial carrier fluxes and the exterior tracer flux are zero.
"""
function sector_grid(nr::Int, na::Int; radius_m, depth_m, omega_rad_s, diffusivity_m2_s)
    nr >= 2 && na >= 4 || throw(ArgumentError("Require at least 2 rings and 4 sectors"))
    all(isfinite, (radius_m, depth_m, omega_rad_s, diffusivity_m2_s)) &&
        radius_m > 0 && depth_m > 0 && diffusivity_m2_s >= 0 ||
        throw(ArgumentError("Require finite geometry and flow, positive sizes, and nonnegative diffusivity"))
    dr = radius_m/nr
    da = 2π/na
    edges = collect(range(0.0, radius_m; length=nr+1))
    inner, outer = edges[1:end-1], edges[2:end]
    centers = (inner+outer)/2
    angles = [(j-.5)*da for j in 1:na]
    radial_area = (outer.^2-inner.^2)/2
    volume = depth_m*da*radial_area
    # Integral of uθ = Ω₀ r (1-r²/R²) over a constant-angle face.
    q = depth_m*omega_rad_s*radial_area .* (1 .- (outer.^2+inner.^2)/(2radius_m^2))
    g_angular = diffusivity_m2_s*depth_m*dr ./ (centers*da)
    g_radial = diffusivity_m2_s*depth_m*da/dr .* edges
    g_radial[end] = 0.0
    (; nr, na, edges, centers, angles, volume, q, g_angular, g_radial)
end

function radial_factor(grid, duration)
    lower = grid.g_radial[2:end-1] ./ grid.volume[2:end]
    upper = grid.g_radial[2:end-1] ./ grid.volume[1:end-1]
    sum_rates = (grid.g_radial[1:end-1]+grid.g_radial[2:end]) ./ grid.volume
    lu(Tridiagonal(-duration/2 .* lower, 1 .+ duration/2 .* sum_rates, -duration/2 .* upper))
end

function angular_factor(grid, i, duration)
    n = grid.na
    diff_rate = grid.g_angular[i]/grid.volume[i]
    adv_rate = grid.q[i]/(2grid.volume[i])
    lower = -duration/2*(diff_rate+adv_rate)
    upper = -duration/2*(diff_rate-adv_rate)
    diagonal = 1+duration*diff_rate
    lhs = spdiagm(-1=>fill(lower,n-1), 0=>fill(diagonal,n), 1=>fill(upper,n-1))
    lhs[1,n] = lower
    lhs[n,1] = upper
    lu(lhs)
end

# Solve each CN update in increment form. Its RHS is the difference of
# shared face fluxes, and vanishes exactly for a spatially uniform field.
function radial_step!(c, rhs, grid, factor, duration)
    @inbounds for j in 1:grid.na
        incoming = 0.0
        for i in 1:grid.nr-1
            outgoing = grid.g_radial[i+1]*(c[i,j]-c[i+1,j])
            rhs[i] = duration*(incoming-outgoing)/grid.volume[i]
            incoming = outgoing
        end
        rhs[end] = duration*incoming/grid.volume[end]
        ldiv!(factor, rhs)
        for i in 1:grid.nr
            c[i,j] += rhs[i]
        end
    end
end

@inline face_flux(left, right, q, g) = q*(left+right)/2+g*(left-right)

function angular_step!(c, rhs, increment, grid, factors, duration)
    @inbounds for i in 1:grid.nr
        q, g = grid.q[i], grid.g_angular[i]
        seam = face_flux(c[i,end], c[i,1], q, g)
        incoming = seam
        scale = duration/grid.volume[i]
        for j in 1:grid.na-1
            outgoing = face_flux(c[i,j], c[i,j+1], q, g)
            rhs[j] = scale*(incoming-outgoing)
            incoming = outgoing
        end
        rhs[end] = scale*(incoming-seam)
        ldiv!(increment, factors[i], rhs)
        for j in 1:grid.na
            c[i,j] += increment[j]
        end
    end
end

function record!(history, c, grid, time)
    lo, hi = extrema(c)
    isfinite(lo) && isfinite(hi) || error("Non-finite tracer concentration")
    mass = sum(sum(c[i,j]*grid.volume[i] for i in 1:grid.nr) for j in 1:grid.na)
    push!(history, (time, mass, lo, hi))
end

"""Conservative full-field transport on the disk at requested physical times."""
function transport(grid, initial, times; dt_max_s)
    size(initial) == (grid.nr,grid.na) || throw(DimensionMismatch("Expected (rings,sectors) field"))
    isfinite(dt_max_s) && dt_max_s > 0 || throw(ArgumentError("Require positive finite time step cap"))
    c = Matrix{Float64}(initial)
    history = NTuple{4,Float64}[]
    record!(history, c, grid, 0.0)
    fields = Matrix{Float64}[]
    previous = 0.0
    radial_rhs = zeros(grid.nr)
    angular_rhs = zeros(grid.na)
    angular_increment = similar(angular_rhs)
    for time in times
        isfinite(time) && time > previous || throw(ArgumentError("Observation times must increase"))
        steps = ceil(Int, (time-previous)/dt_max_s)
        dt = (time-previous)/steps
        radial = radial_factor(grid, dt/2)
        angular = [angular_factor(grid, i, dt) for i in 1:grid.nr]
        for step in 1:steps
            radial_step!(c, radial_rhs, grid, radial, dt/2)
            angular_step!(c, angular_rhs, angular_increment, grid, angular, dt)
            radial_step!(c, radial_rhs, grid, radial, dt/2)
            record!(history, c, grid, step == steps ? time : previous+step*dt)
        end
        push!(fields, copy(c))
        previous = time
    end
    (; fields, history)
end

end
