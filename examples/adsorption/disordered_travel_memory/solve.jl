include("model.jl")
using .TravelMemory: predict
using DelimitedFiles

length(ARGS) == 3 || error("Use training.f64 probe_passages.f64 output_directory")
Base.ENDIAN_BOM == 0x04030201 || error("Expected little-endian storage")
load_arrivals(path,columns) = transpose(reshape(reinterpret(Float64,read(path)),columns,:))
train, probes = load_arrivals(ARGS[1],2), load_arrivals(ARGS[2],8)
out = ARGS[3]
n = size(train,1)
size(probes,1) == 2n || error("Expected two independent probe halves")
for (name, groups, passages) in (("baseline",1,probes),("two",2,probes),("four",4,probes),
        ("eight",8,probes),("four_first",4,view(probes,1:n,:)),("four_second",4,view(probes,n+1:2n,:)))
    model = predict(train,passages,groups)
    writedlm(joinpath(out,"$name.csv"),model.prediction,',')
    if name == "four"
        open(joinpath(out,"rules.f64"),"w") do file
            for rule in model.rules
                write(file,rule.transition,rule.delay,rule.delay2)
            end
        end
    end
end
