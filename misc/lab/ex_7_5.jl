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

function idft(X)
    
    N = length(X)
    x = zeros(ComplexF64,N)

    for k in 0:N-1
        for n in 0:N-1
            x[n+1] += X[k+1] * exp(im * 2pi * k * n / N) 
        end
    end
    return (1/N) .* x
end

signal = sin.(-2pi:0.01:2pi)

sin_dft = dft(signal)
sin_idft = idft(sin_dft)

plot(real.(sin_idft))