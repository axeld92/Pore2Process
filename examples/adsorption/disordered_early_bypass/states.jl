module BypassStates

include("../disordered_travel_memory/model.jl")
using .TravelMemory: age_sizes, states_by_age
using .TravelMemory.SectionExchange: training, paths

function classify(train, pass)
    fit = training(train)
    sizes = age_sizes(fit.sizes,4)
    scale, remainder = divrem(size(pass,1),size(train,1))
    scale>0 && remainder==0 && size(pass,2)==8 || error("Invalid tracer archive split")
    delays = diff(hcat(zeros(size(pass,1)),pass);dims=2)
    hcat([states_by_age(paths(delays[:,k],fit.sizes),pass[:,k],sizes,scale) for k in 1:8]...)
end

end
