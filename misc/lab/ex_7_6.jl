using CairoMakie

function goertzel(x,k)
    N = length(x)
    s_p = 0
    s_pp = 0
        for n in 0:N-1
            s = x[n+1] + 2 * cos(2pi/N) * s_p - s_pp
            s_pp = s_p
            s_p = s
        end
        X = s_p - exp(-im * 2pi * k / N) * s_pp
    return X
end

signal = sin.(-2pi:0.01:2pi)

X = goertzel(signal, 5)