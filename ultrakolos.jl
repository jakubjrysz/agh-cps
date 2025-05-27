using CairoMakie
using LinearAlgebra

mean(x::AbstractVector)::Real = sum(x)/length(x)

energy(x::AbstractVector)::Real = sum(abs2, x)

power(x::AbstractVector)::Real = energy(x)/length(x)

rms(x::AbstractVector)::Real = sqrt(power(x))




ramp_wave(t::Real)::Real = 2 * (t % 1)

triangle_wave(t::Real)::Real = 4*(0.5 - abs((t % 1) - 0.5))-1

sawtooth_wave(t::Real)::Real = 2 * (t % 1) - 1

square_wave(t::Real)::Real = ifelse(mod(t,1) < 0.5, 1 , -1)


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

function lti_amp(b::Vector, a::Vector, F::Vector)
                                                    #z = exp(j2π*f)     gdzie f = f/fp 
    f_len = length(F)
    b_len = length(b)
    a_len = length(a)
    licznik = zeros(ComplexF64, f_len)
    mianownik = ones(ComplexF64, f_len)

    for i in 1:f_len
        z = exp(im * 2 * π * F[i])
        for j in 1:b_len
            licznik[i] += b[j] * z ^ (- j + 1)
        end
        for k in 2:a_len
            mianownik[i] += a[k] * z ^ (- k + 1)
        end
        amp[i] = abs(licznik[i]/mianownik[i])
    end
    return amp
end

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

function interpolate(s::Vector, m::Vector, t::Vector)
    delta_t = m[2] - m[1]
    x = [sum(s[j] * sinc((ti - m[j]) / delta_t) for j in eachindex(m))  for ti in t]
    return x
end

function quantize(x, L)
    y = L[argmin(abs.(x .- L))]
    return y
end

function quantization_error(a, b, x)
    N = 5   #liczba bitów
    L = range(start = a, stop = b, length = 2 ^ N)
    x_quantized = [quantize(xi, L) for xi in x]
    error = x .- x_quantized
    return error
end