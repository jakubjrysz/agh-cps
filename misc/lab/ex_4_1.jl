using Random
using CairoMakie
using LinearAlgebra


function triangular_wave(t)
    return 2*abs(2*((4*t+0.25) - floor(0.75+4*t)))-1
end

function mean(x)
    l = length(x)
    s = sum(x)
    fac = 1/(2*l+1)
    return fac*sum(x)
end


x = 0:0.001:0.1

lines(triangular_wave.(x))
mean(triangular_wave.(x))