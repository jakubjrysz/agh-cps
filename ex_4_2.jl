using Random
using CairoMakie
using LinearAlgebra


function peak2peak(x)
    return maximum(x) - minimum(x)
end


x = 0:0.001:10

lines(sin.(x))
peak2peak(sin.(x))