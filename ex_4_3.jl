using Random
using CairoMakie
using LinearAlgebra


function energy(x)
    return sum(x.^2)
end


x = 0:0.001:1

lines(sin.(x))
energy(sin.(x))