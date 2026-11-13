module DisorderedExchange

using LinearAlgebra

function deciles(durations)
    n = length(durations)
    limits = round.(Int,n .* (0:10) ./ 10)
    labels = Vector{Int}(undef,n)
    order = sortperm(durations)
    for c in 1:10
        labels[order[limits[c]+1:limits[c+1]]] .= c
    end
    labels
end

function calibrate(arrivals,hits)
    n = size(arrivals,1)
    size(arrivals)==size(hits)==(n,2) && n>=10 ||
        throw(ArgumentError("Expected the first two sections and at least ten markers"))
    duration = arrivals[:,2]-arrivals[:,1]
    all(isfinite,arrivals) && all(>(0),arrivals[:,1]) && all(>(0),duration) ||
        throw(ArgumentError("Training arrivals must be finite, complete and increasing"))
    all((isnan(h) || (0<=h<=arrivals[i,1])) for (i,h) in enumerate(hits[:,1])) &&
        all((isnan(h) || (arrivals[i,1]<=h<=arrivals[i,2])) for (i,h) in enumerate(hits[:,2])) ||
        throw(ArgumentError("Training wall clocks fall outside their sections"))
    first,second = deciles(arrivals[:,1]),deciles(duration)
    counts = zeros(Int,10,10)
    clear = zeros(Int,10,10)
    delay = zeros(10,10)
    delay2 = zeros(10,10)
    states = zeros(10,6)
    histories = zeros(2,2)
    for k in 1:n
        i,j = first[k],second[k]
        counts[i,j] += 1
        clear[i,j] += isnan(hits[k,2])
        delay[i,j] += duration[k]
        delay2[i,j] += duration[k]^2
        t = arrivals[k,2]
        alive1 = isnan(hits[k,1])
        alive2 = alive1 && isnan(hits[k,2])
        fast = i<=9
        states[j,1] += 1
        states[j,2] += t
        states[j,3] += t^2
        states[j,4] += alive2
        states[j,5] += alive2 && fast
        states[j,6] += fast && j<=9
        histories[1,1] += alive1
        histories[2,1] += alive2
        histories[1,2] += alive1 && fast
        histories[2,2] += alive2 && fast
    end
    donor = vec(sum(counts;dims=2))
    receiver = vec(sum(counts;dims=1))
    donor==receiver || error("Carrier class fractions do not match")
    totals = [n,sum(donor[1:9])]
    transition = counts./donor
    wallfree = clear./donor
    receiver_survival = vec(sum(clear;dims=1))./receiver
    factorized = transition.*transpose(receiver_survival)
    flux = donor/n
    maximum(abs,transpose(transition)*flux-flux)<=1e-14 ||
        error("Uniform inert concentration is not preserved")
    (;counts,clear,delay=delay./max.(counts,1),delay2=delay2./max.(counts,1),
      transition,wallfree,factorized,flux,states=states/n,
      histories=histories./transpose(totals),fast_share=totals[2]/n)
end

"""Advance inert moments and wall-free dose with frozen donor-by-receiver rules."""
function propagate(model,periods)
    mass,moment,second,alive,fast,persistent = eachcol(copy(model.states))
    # Holdout rows contain total and fast survival, their conditional capture,
    # inert mean and variance, and uninterrupted fast-path persistence.
    results = zeros(periods-2,7)
    survival = sum(alive)
    fast_survival = sum(fast)/model.fast_share
    T = transpose(model.transition)
    W = transpose(model.wallfree)
    D = transpose(model.transition.*model.delay)
    D2 = transpose(model.transition.*model.delay2)
    for row in axes(results,1)
        second = T*second+2D*moment+D2*mass
        moment = T*moment+D*mass
        mass = T*mass
        alive = W*alive
        fast = W*fast
        persistent = T*persistent
        persistent[end] = 0.
        a,b = sum(alive),sum(fast)/model.fast_share
        mean = sum(moment)
        isapprox(sum(mass),1.;atol=1e-12,rtol=0.) &&
            0<=a<=survival+1e-12 && 0<=b<=fast_survival+1e-12 &&
            sum(second)>=mean^2 || error("Tracer balance or bounds failed")
        results[row,:] = [a, b, 1-a/survival, 1-b/fast_survival, mean,
                          sum(second)-mean^2, sum(persistent)/model.fast_share]
        survival,fast_survival = a,b
    end
    results
end

end
