using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function isPrime(x)
    for i in 2:x-1
        if x%i == 0
            return "0"
            break
        else
            continue
        end
    end
    return "1"
end

isPrime(6)