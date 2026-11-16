module ColumnHandoff

include("../disordered_travel_memory/model.jl")
using .TravelMemory: age_sizes, states_by_age, memory_rule, predict
using .TravelMemory.SectionExchange: training, paths, rule
using ..MemoryColumn: delays

function load_transport(training_path,probe_path,rule_path,groups)
    Base.ENDIAN_BOM == 0x04030201 || error("Expected little-endian storage")
    load(path,columns) = transpose(reshape(reinterpret(Float64,read(path)),columns,:))
    train,raw = load(training_path,2),load(probe_path,16)
    passage = raw[:,1:8]
    fit = training(train)
    sizes = age_sizes(fit.sizes,groups)
    flow = vec(sizes)/sum(sizes)
    n = length(flow)
    probesize, remainder = divrem(size(raw,1),size(train,1))
    probesize==2 && remainder==0 || error("Unexpected probe cohort size")
    late = predict(train,passage,groups).rules
    if groups==4
        saved = reshape(reinterpret(Float64,read(rule_path)),n,n,3,6)
        for (k,r) in enumerate(late)
            actual = (r.transition,r.delay,r.delay2)
            for (j,table) in enumerate(actual)
                table == saved[:,:,j,k] || error("Age-conditioned rule differs from the qualified run")
            end
        end
        late = [(transition=saved[:,:,1,k],delay=saved[:,:,2,k],delay2=saved[:,:,3,k]) for k in 1:6]
    end
    age1,age2 = train[:,1],train[:,2]
    duration2 = age2-age1
    labels(d,a,scale) = states_by_age(paths(d,fit.sizes),a,sizes,scale)
    first,second = labels(age1,age1,1),labels(duration2,age2,1)
    delay1,delay1sq = zeros(n),zeros(n)
    for i in eachindex(age1)
        delay1[first[i]] += age1[i]/length(age1)
        delay1sq[first[i]] += age1[i]^2/length(age1)
    end
    start = delays(reshape(flow,1,:),reshape(delay1,1,:),reshape(delay1sq,1,:),[1.])
    next = memory_rule(first,second,duration2,sizes,fit.stationary)
    sections = [start,[delays(r.transition,r.delay,r.delay2,flow) for r in [next;late]]...]
    # Separate contact data from the independent operator probes.
    local_times = diff(hcat(zeros(size(raw,1)),passage);dims=2)
    hits = zeros(Int,n,8)
    for k in 1:8
        state = labels(local_times[:,k],passage[:,k],probesize)
        for i in eachindex(state)
            event = raw[i,k+8]
            if !isnan(event)
                lo = k==1 ? 0. : passage[i,k-1]
                lo <= event <= passage[i,k] || error("Wall hit falls outside its section")
                hits[state[i],k] += 1
            end
        end
    end
    totals = sum(hits;dims=1)
    minimum(totals)>0 || error("Every section needs sampled solid-wall encounters")
    area = hits./totals
    tau = [sum(r.share.*r.seconds) for r in sections]
    (; sections, flow, area, hits, tau)
end

end
