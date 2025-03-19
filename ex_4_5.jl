using Random
using CairoMakie
using LinearAlgebra


function rms(x)
    n = length(x)
    return sqrt((sum(x.^2))/n)
end


x = -10*pi:0.001:10*pi

y = sin.(x)
lines(sin.(x))
println(sin.(x))

