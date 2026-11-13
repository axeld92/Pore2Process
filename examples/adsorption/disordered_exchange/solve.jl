include("model.jl")
using .DisorderedExchange: calibrate, propagate
using DelimitedFiles

length(ARGS)==2 || error("Use training.f64 output_directory")
train = reshape(reinterpret(Float64,read(ARGS[1])),4,:)
model = calibrate(transpose(train[1:2,:]),transpose(train[3:4,:]))
output = ARGS[2]
for key in (:counts,:clear,:delay,:delay2,:transition,:wallfree,:factorized,:flux,:states,:histories,:fast_share)
    value = getproperty(model,key)
    writedlm(joinpath(output,"$key.csv"),value isa Number ? [value] : value,',')
end
writedlm(joinpath(output,"pair.csv"),propagate(model,8),',')
writedlm(joinpath(output,"receiver.csv"),propagate(merge(model,(wallfree=model.factorized,)),8),',')
