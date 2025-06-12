using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function ci_rectangular(T,t)
    if t<T/2 && t>-T/2
        return 1
    elseif t == T/2 || t == -T/2
        return 0.5
    else
        return 0
    end
end

ci_rectangular(10,5)
