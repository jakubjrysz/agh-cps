using CairoMakie
function dft(x)
    
    N = length(x)
    X = zeros(ComplexF64,N)

    for k in 0:N-1
        for n in 0:N-1
            X[k+1] += x[n+1] * exp(-im * 2pi * k * n / N) 
        end
    end
    return X
end

signal = sin.(-2pi:0.01:2pi)

X = dft(signal)

plot(abs.(X))