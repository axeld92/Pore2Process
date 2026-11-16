module TravelMemory

include("../disordered_section_exchange/model.jl")
using .SectionExchange: training, paths, local_rules, propagate

function age_sizes(sizes, groups)
    groups > 0 && groups <= minimum(sizes) || throw(ArgumentError("Invalid age group count"))
    hcat([diff(div.(n .* (0:groups),groups)) for n in sizes]...)
end

"""Index the group in age order within its local path."""
function states_by_age(path, age, sizes, scale)
    groups, classes = size(sizes)
    labels = similar(path)
    for i in 1:classes
        members = findall(==(i),path)
        edges = [0; cumsum(scale .* sizes[:,i])]
        length(members) == last(edges) || error("Unexpected class population")
        order = sortperm(age[members])
        for g in 1:groups
            labels[members[order[edges[g]+1:edges[g+1]]]] .= (i-1)*groups+g
        end
    end
    labels
end

function initial_states(arrival, fit, sizes)
    age = arrival[:,2]
    labels = states_by_age(paths(age-arrival[:,1],fit.sizes),age,sizes,1)
    initial = zeros(length(sizes),3)
    for k in eachindex(age)
        initial[labels[k],:] .+= (1.,age[k],age[k]^2)
    end
    initial ./= length(age)
    groups, classes = size(sizes)
    collapsed = dropdims(sum(reshape(initial,groups,classes,3);dims=1);dims=1)
    isapprox(collapsed,fit.states;rtol=1e-13,atol=1e-15) ||
        error("Age grouping changed the section 2 calibration")
    initial
end

function memory_rule(donor, receiver, delay, sizes, baseline)
    groups, classes = size(sizes)
    counts = zeros(Int,length(sizes),length(sizes))
    first = zeros(size(counts))
    second = zeros(size(counts))
    for k in eachindex(delay)
        i,j = donor[k],receiver[k]
        counts[i,j] += 1
        first[i,j] += delay[k]
        second[i,j] += delay[k]^2
    end
    rows = vec(sum(counts;dims=2))
    columns = vec(sum(counts;dims=1))
    fractions = vec(sizes)/sum(sizes)
    rows == columns && isapprox(rows/sum(rows),fractions;rtol=1e-14,atol=0.) ||
        error("Age-dependent exchange changed steady state fractions")

    collapse(values) = dropdims(sum(reshape(values,groups,classes,groups,classes);dims=(1,3));dims=(1,3))
    path_rows = vec(sum(baseline.counts;dims=2))
    collapse(counts) == baseline.counts &&
        isapprox(collapse(first)./path_rows,baseline.delay;rtol=1e-13,atol=0.) &&
        isapprox(collapse(second)./path_rows,baseline.delay2;rtol=1e-13,atol=0.) ||
        error("Age grouping changed the section-specific exchange or local delay")
    (;counts,transition=counts./rows,delay=first./rows,delay2=second./rows)
end

function predict(arrival, passages, groups)
    fit = training(arrival)
    sizes = age_sizes(fit.sizes,groups)
    size(passages,2) == 8 || throw(ArgumentError("Expected probe passages for sections 1–8"))
    scale, remainder = divrem(size(passages,1),size(arrival,1))
    scale > 0 && remainder == 0 || throw(ArgumentError("Expected a multiple of the training count"))
    all(x -> isfinite(x) && x > 0, passages[:,1]) ||
        throw(ArgumentError("Expected complete, positive first probe passages"))
    durations = diff(passages;dims=2)
    baseline = local_rules(durations,fit.sizes)
    baseline_prediction = propagate(fit.states,baseline)
    labels = [states_by_age(paths(d,fit.sizes),a,sizes,scale)
              for (d,a) in zip(eachcol(durations),eachcol(view(passages,:,2:8)))]
    rules = [memory_rule(labels[k],labels[k+1],durations[:,k+1],sizes,baseline[k])
             for k in eachindex(baseline)]
    prediction = propagate(initial_states(arrival,fit,sizes),rules)
    groups == 1 && !isapprox(prediction,baseline_prediction;rtol=1e-13,atol=1e-15) &&
        error("One age group did not reproduce the baseline")
    (; prediction, baseline_prediction, rules)
end

end
