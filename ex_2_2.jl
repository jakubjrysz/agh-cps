using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function fact(x)
    y = 1
    for i in 2:1:x
        y *= i
    end
    return y
end

fact(6)