using Random
using CairoMakie
using LinearAlgebra


function power(x)
    fac = 1/(lastindex(x) - firstindex(x)+1)
    return fac*sum(x.^2)
end


x = -2*pi:0.001:2*pi

lines(sin.(x))
power(sin.(x))

lastindex(x) - firstindex(x) + 1
println(length(x))