module RepeatJoin

using LinearAlgebra
using ..MemoryColumn: Delays, solve_column

function path_preserving(rules, blocks)
    blocks >= 1 || throw(ArgumentError("Expected at least one repeat"))
    first = rules[1]
    all(==(1),first.donor) || throw(ArgumentError("The baseline first section must use one feed"))
    # Each old first-section branch retains its receiving state, flux and delay.
    join = Delays(first.receiver,first.receiver,first.share,first.seconds)
    [rules;repeat([join;rules[2:end]],blocks-1)]
end

function solve_preserved(rules,flow,contact;blocks,kwargs...)
    solve_column(path_preserving(rules,blocks),flow,repeat(contact,1,blocks);
                 blocks=1,kwargs...)
end

"""First two moments of accumulated delay, from the discrete network."""
function delay_moments(rules,flow)
    n = length(flow)
    mass, first, second = zeros(n),zeros(n),zeros(n)
    mass[1] = 1.
    for rule in rules
        m,f,s = zeros(n),zeros(n),zeros(n)
        donor_flow = zeros(n)
        for a in eachindex(rule.share)
            donor_flow[rule.donor[a]] += rule.share[a]
        end
        for a in eachindex(rule.share)
            i,j,t = rule.donor[a],rule.receiver[a],rule.seconds[a]
            p = rule.share[a]/donor_flow[i]
            m[j] += p*mass[i]
            f[j] += p*(first[i]+t*mass[i])
            s[j] += p*(second[i]+2t*first[i]+t^2*mass[i])
        end
        isapprox(m,flow;rtol=1e-12,atol=1e-15) || error("Water flow changed at a join")
        mass,first,second = m,f,s
    end
    mean = sum(first)
    (;mean,variance=sum(second)-mean^2)
end

function repeat_moments(rules,flow,blocks,holding_s,preserved)
    block = delay_moments(rules,flow)
    scale = holding_s/(blocks*block.mean)
    if preserved
        m = delay_moments(path_preserving(rules,blocks),flow)
        return (;mean=m.mean*scale,variance=m.variance*scale^2)
    end
    (;mean=blocks*block.mean*scale,variance=blocks*block.variance*scale^2)
end

end
