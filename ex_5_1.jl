using Random
using CairoMakie
using LinearAlgebra


function quantize(x,l)
    
    q = zeros(length(x))
    trueQ = zeros(length(x))

    for i in eachindex(x)
        q[i] = argmin((abs.(l .- x[i])))
        trueQ[i] = l[convert(Int,q[i])] 
    end
    return trueQ 
end

x = -2*pi:0.1:2*pi
l = -1:0.5:1

t1 = time()

quantize(sin.(x),l)

elapsed_time = time() - t1

println(elapsed_time)
scatter(quantize(sin.(x),l))

