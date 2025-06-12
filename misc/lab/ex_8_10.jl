using CairoMakie
using FFTW

function amplitude_spectrum(x, w)

    x_windowed = w .* x
    X = fft(x_windowed)
    A = abs.(X)

    return A

end

t = 0:0.01:10

x = 2*sin.(30t) + cos.(25t)

amplitude_spectrum(x,0.8)
plot(amplitude_spectrum(x,0.8))