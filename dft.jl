using CairoMakie
using LinearAlgebra
using FFTW

function myDFT(signal)
    n = length(signal)
    result = zeros(Complex{Float64}, n)
    for k in 0:n-1
        for j in 0:n-1
            result[k + 1] += signal[j + 1] * exp(-2im * π * k * j / n)
        end
    end
    return result
end

signal = [0.85 - 0.19im, 0.12 + 0.38im, -0.39 - 0.74im, -0.43 - 0.74im, -0.93 - 0.38im, -1.9 + 0.58im, 0.34 - 0.62im, 0.21 + 0.22im, 0.5 + 1.05im, -0.15 + 1.26im, -0.7 + 0.02im, -0.89 - 0.26im, -0.1 - 0.09im, 0.78 + 0.49im, 0.05 + 1.0im, -1.49 - 1.07im, 0.05 - 1.5im, 0.36 + 0.28im, -0.32 + 1.51im, 0.61 + 0.15im, 0.09 - 0.07im, 0.81 - 0.27im, 0.05 + 0.09im, 0.17 + 0.83im, -0.59 + 0.36im, -0.8 - 0.15im, -0.57 + 0.21im, -0.2 - 0.84im, -0.3 - 0.13im, -0.79 + 0.12im, -0.79 + 0.65im, 0.7 - 0.66im, 0.25 + 0.6im, 0.26 - 0.71im, -0.36 + 1.58im, 0.56 - 0.1im, -0.45 - 0.23im, 0.03 + 0.51im, -0.76 - 1.61im, 1.63 + 0.98im, -0.83 - 0.48im]
signal_dft = myDFT(signal)

fp = 861
f_res = fp/length(signal)

lines(abs.(signal_dft))