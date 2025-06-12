using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function isEven(x)

    x%2 == 0 ? "1" : "0"

end

isEven(6)