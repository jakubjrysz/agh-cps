using Random
using CairoMakie
using LinearAlgebra


function power(x)
    fac = 1/(length(x))
    return fac*sum(x.^2)
end


x = -2*pi:0.001:2*pi

power(sin.(x))