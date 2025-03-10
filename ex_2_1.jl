using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function fact(x)
    if x == 0
        return 1
    end

    return x*fact(x-1)
end

fact(5)