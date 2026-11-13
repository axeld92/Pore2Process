module SectionExchange

include("../disordered_exchange/model.jl")
using .DisorderedExchange: deciles

const CLASSES = 10

"""Rank-based paths with the calibrated class sizes, scaled to the probe count."""
function paths(duration, sizes)
    scale, remainder = divrem(length(duration), sum(sizes))
    scale > 0 && remainder == 0 || error("Probe count must be a multiple of the training count")
    edges = cumsum([0; scale .* sizes])
    order = sortperm(duration)
    labels = similar(order)
    for i in eachindex(sizes)
        labels[order[edges[i]+1:edges[i+1]]] .= i
    end
    labels
end

function rule(donor, receiver, delay, sizes)
    counts = zeros(Int, CLASSES, CLASSES)
    first = zeros(CLASSES, CLASSES)
    second = zeros(CLASSES, CLASSES)
    for k in eachindex(delay)
        i, j = donor[k], receiver[k]
        counts[i,j] += 1
        first[i,j] += delay[k]
        second[i,j] += delay[k]^2
    end
    total = vec(sum(counts; dims=2))
    flux = sizes/sum(sizes)
    transition = counts ./ total
    maximum(abs, transpose(transition)*flux - flux) <= 1e-14 ||
        error("Transition changed the steady class flux fractions")
    # Probability-weighted delay moments. Unsampled pairs carry zero probability.
    (; counts, transition, delay=first ./ total, delay2=second ./ total)
end

function training(arrival)
    size(arrival,2) == 2 && size(arrival,1) >= CLASSES ||
        throw(ArgumentError("Use only sections 1–2 for training"))
    n = size(arrival,1)
    t = arrival[:,2]
    d1, d2 = arrival[:,1], t - arrival[:,1]
    all(x -> isfinite(x) && x > 0, d1) && all(x -> isfinite(x) && x > 0, d2) ||
        throw(ArgumentError("Training passages must be complete and strictly increasing"))
    a, b = deciles(d1), deciles(d2)
    sizes = [count(==(i), b) for i in 1:CLASSES]
    states = zeros(CLASSES,3)
    for k in eachindex(t)
        states[b[k],:] .+= (1., t[k], t[k]^2)
    end
    (; states=states/n, sizes, stationary=rule(a,b,d2,sizes))
end

function local_rules(durations, sizes)
    size(durations,2) == 7 || throw(ArgumentError("Use probe durations for sections 2–8"))
    all(x -> isfinite(x) && x > 0, durations) ||
        throw(ArgumentError("Probe passages must be complete and strictly increasing"))
    labels = [paths(d,sizes) for d in eachcol(durations)]
    [rule(labels[k-1],labels[k],durations[:,k],sizes) for k in 2:length(labels)]
end

function propagate(initial, rules)
    mass, age, age2 = eachcol(initial)
    output = zeros(length(rules),3)
    for (k, r) in enumerate(rules)
        T, D, D2 = transpose(r.transition), transpose(r.delay), transpose(r.delay2)
        age2 = T*age2 + 2D*age + D2*mass
        age = T*age + D*mass
        mass = T*mass
        mean = sum(age)
        variance = sum(age2) - mean^2
        isapprox(sum(mass), 1.; atol=1e-12, rtol=0.) && all(>=(0),mass) &&
            isfinite(variance) && variance >= 0 || error("Tracer balance or bounds failed")
        output[k,:] .= (sum(mass),mean,variance)
    end
    output
end

end
