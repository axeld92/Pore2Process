include("model.jl")
using .SectionExchange: training, local_rules, propagate
using DelimitedFiles

length(ARGS) == 3 || error("Use training.f64 probe_durations.f64 output_directory")
Base.ENDIAN_BOM == 0x04030201 || error("Expected little-endian storage")
arrival = transpose(reshape(reinterpret(Float64,read(ARGS[1])),2,:))
duration = transpose(reshape(reinterpret(Float64,read(ARGS[2])),7,:))
output = ARGS[3]
fit = training(arrival)
n = size(arrival,1)
size(duration,1) == 2n || error("Expected two independent probe sets of training size")
sections = local_rules(duration,fit.sizes)
writedlm(joinpath(output,"stationary.csv"),propagate(fit.states,fill(fit.stationary,6)),',')
writedlm(joinpath(output,"sectional.csv"),propagate(fit.states,sections),',')
for (name, rows) in (("first",1:n),("second",n+1:2n))
    rules = local_rules(view(duration,rows,:),fit.sizes)
    writedlm(joinpath(output,"probe_$name.csv"),propagate(fit.states,rules),',')
end
for (i, section) in enumerate(sections), name in (:counts,:transition,:delay,:delay2)
    writedlm(joinpath(output,"section_$(i+2)_$name.csv"),getproperty(section,name),',')
end
