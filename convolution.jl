using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

function conv(signal,kernel)
    m = length(signal)
    n = length(kernel)
    reversedKernel = reverse(kernel)
    outputLength = m + n - 1
    result = zeros(outputLength)

    for i in 1:outputLength
        for j in 1:n
            if i - j + 1 ≥ 1 && i - j + 1 ≤ m
                result[i] += signal[i - j + 1] * reversedKernel[j]  
            end
        end
    end
    return result
end

x::Vector{Float64} = [-0.16, 1.83, -3.23, 2.58, 2.77, 0.74, 3.55, -4.19, -2.08, 1.2, 3.27, -2.16, 1.49, -0.4, -0.42, -2.84, -2.42, -4.6, 0.98, 3.51, 3.42, 2.81, 3.97, 1.1, 2.5, 1.33, 1.22, -3.13, 0.73, -0.27, -2.01, -4.18, -2.38, 1.74, -3.92, -2.47, -0.96, 0.59, 4.34, -1.02, 2.31, -1.81, -4.0, 2.98, 0.57, 0.7, -0.97, 4.5, 1.97, 0.76, -1.52, 4.36, 4.69, -2.78, -0.71, 1.0, 1.48, -3.53, -0.87, -2.57, -3.13, -1.59, 2.74]
h::Vector{Float64} = [2.13, 1.32, 3.18, -0.43, 0.68, 3.58, -0.06, -0.67, 0.83, -1.1, -0.97, 1.92, -2.09, -0.09, -3.91]

y = conv(x,h)

pwr = sum(y.^2)/length(y)