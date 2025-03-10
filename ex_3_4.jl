using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function ci_triangular(T,t)
    if t<T/2 && t>-T/2
        return (0.5*T-abs(t))/(0.5*T)
    else
        return 0
    end
end

ci_triangular(10,2)
