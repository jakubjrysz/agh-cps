module CPS

using CairoMakie
using LinearAlgebra
using OffsetArrays
using FFTW
using SpecialFunctions


const author = Dict(
    "index" => "421266",
    "name" => "Jakub Rysz",
    "email" => "jrysz@student.agh.edu.pl",
)



###############################################################################
# Parametry sygnałów                                                          #
###############################################################################
mean(x::AbstractVector)::Number = sum(x)/length(x)
peak2peak(x::AbstractVector)::Real = abs(maximum(x) - minimum(x))
energy(x::AbstractVector)::Real = sum(x.^2)
power(x::AbstractVector)::Real = sum(x.^2)/length(x)
rms(x::AbstractVector)::Real = sqrt(power(x))

function running_mean(x::AbstractVector, M::Integer)::Vector
    result::AbstractVector = zeros(length(x))
    for k in 1:length(x)
        n1 = k - M < 1 ? 1 : k - M
        n2 = k + M > lastindex(x) ? lastindex(x) : k + M
        result[k] = (1 / (n2 - n1 + 1)) * sum(x[n1:n2])
    end
    return result
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    result::AbstractVector = zeros(length(x))
    for k in 1:length(x)
        n1 = k - M < 1 ? 1 : k - M
        n2 = k + M > lastindex(x) ? lastindex(x) : k + M
        result[k] = sum(abs2, x[n1:n2])
    end
    return result
end

function running_power(x::AbstractVector, M::Integer)::Vector
    result::AbstractVector = zeros(length(x))
    for k in 1:length(x)
        n1 = k - M < 1 ? 1 : k - M
        n2 = k + M > lastindex(x) ? lastindex(x) : k + M
        result[k] = (1 / (n2 - n1 + 1)) * sum(abs2, x[n1:n2])
    end
    return result
end



###############################################################################
# Modelowanie sygnałów                                                        #
###############################################################################
cw_rectangular(t::Real; T=1.0)::Real = abs(t) < T / 2 ? 1 : (abs(t) == T / 2 ? 0.5 : 0)
cw_triangle(t::Real; T=1.0)::Real = abs(t) < T ? 1 - abs(t) : 0
cw_literka_M(t::Real; T=1.0)::Real = abs(t) < T ? (t < 0 ? -t + 1 : t + 1) : 0
cw_literka_U(t::Real; T=1.0)::Real = abs(t) < T ? t^2 : 0

ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
triangular_wave(t::Real)::Real = 2 / π * asin(sin(2π * t))
square_wave(t::Real)::Real = ifelse(mod(t, 1) < 0.5, 1, -1)
pulse_wave(t::Real, ρ::Real)::Real = ifelse(mod(t, 1) < ρ, 1, 0)

function impuse_repeater(g::Function, t0::Real, t1::Real)::Function 
    x -> g(mod(x - t1, t2 - t1) + t1)
end

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    temp = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band * 2π
        temp += (-1)^n * sin.(arg * t) / n
        n += 1
    end
    signal += -2A / π * temp
    return signal
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band * 2π
        signal += (-1)^n * sin.(arg * t) / n
        n += 1
    end
    signal *= 2A / π
    return signal
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band * 2π
        signal += (-1)^((n - 1) / 2) * sin.(arg * t) / n^2
        n += 2
    end
    signal *= (8A / π^2)
    return signal
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * (2n - 1) * (1 / T)) < band * 2π
        signal += sin.(arg * t) / (2n - 1)
        n += 1
    end
    signal *= 4 * A / π
    return signal
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    signal = (sawtooth_wave_bl.(t .- (T / 2); A, T, band) - sawtooth_wave_bl.(t .- ((T / 2) + ρ); A, T, band)) .+ (2 * A * ρ)
    return signal
end

function impulse_repeater_bl(g::Function, t1::Real, t2::Real, band::Real)::Function
    T::Float64 = t2 - t1
    ω₀::Float64 = (2π / T)
    n_terms::Integer = div(band * 2π, ω₀)

    N = 1000
    ts = range(t1, t2, length=N + 1)
    Δt = (t2 - t1) / N

    function trapezoidal_integral(f::Function, ts::AbstractVector, Δt::Float64)
        integral = 0.5 * (f(ts[1]) + f(ts[end]))
        for i in 2:N
            integral += f(ts[i])
        end
        return integral * Δt
    end

    a0 = 1 / T * trapezoidal_integral(g, ts, Δt)
    an_coeffs = zeros(Float64, n_terms)
    bn_coeffs = zeros(Float64, n_terms)

    for n in 1:n_terms
        an_coeffs[n] = 2 / T * trapezoidal_integral(t -> g(t) * cos(ω₀ * n * t), ts, Δt)
        bn_coeffs[n] = 2 / T * trapezoidal_integral(t -> g(t) * sin(ω₀ * n * t), ts, Δt)
    end

    function fourier_series_output(t::Float64)
        result = a0 / 2
        for n in 1:n_terms
            result += an_coeffs[n] * cos(ω₀ * n * t) + bn_coeffs[n] * sin(ω₀ * n * t)
        end
        return result
    end

    return fourier_series_output
end

function rand_signal_bl(f1::Real, f2::Real)::Function
    f = f1 .+ rand(1000) .* (f2 - f1)
    ϕ = -π .+ rand(1000) * 2π
    A = randn(1000)
    A = A ./ sqrt(0.5 * sum(A .^ 2))
    return t -> sum(A .* sin.(2π * f .* t .+ ϕ))
end



kronecker(n::Integer)::Real = ifelse(n == 0, 1, 0)
heaviside(n::Integer)::Real = ifelse(n < 0, 0, 1)

# Dyskretne okna czasowe
rect(N::Integer)::AbstractVector{<:Real} = ones(N)
triang(N::Integer)::AbstractVector{<:Real} = [1 - (2abs(n - ((N - 1) / 2))) / (N - 1) for n = 0:N-1]
hanning(N::Integer)::AbstractVector{<:Real} = [0.5(1 - cos(2π * n / (N - 1))) for n = 0:N-1]
hamming(N::Integer)::AbstractVector{<:Real} = [0.54 - 0.46cos(2π * n / (N - 1)) for n = 0:N-1]
blackman(N::Integer)::AbstractVector{<:Real} = [0.42 - 0.5cos(2π * n / (N - 1)) + 0.08cos(4π * n / (N - 1)) for n = 0:N-1]

function chebwin(N::Int, α::Real)

    tg = 10^(α / 20)
    M = N - 1
    n = 0:M
    x = cos.(π * n / M)

    T = cos.(M * acos.(x))

    T ./= maximum(T)

    w = real(ifft(T))
    w = fftshift(w)

    return w[1:N]
end

function kaiser(N::Int, β::Real)
    n = 0:N-1
    α = (N - 1) / 2
    denom = besseli(0, β)

    w = besseli.(0, β * sqrt.(1 .- ((n .- α) ./ α).^2)) ./ denom
    return w
end


###############################################################################
# Próbkowanie i kwantyzacja                                                   #
###############################################################################

quantize(L::AbstractVector)::Function = x -> L[argmin(abs.(x .- L))]
SQNR(N::Integer)::Real =  1.76 + 6.02 * N
SNR(Psignal::Real, Pnoise::Real)::Real = 10 * log10(Psignal / Pnoise)


function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function=sinc
)::Function
    
end


###############################################################################
# Obliczanie dyskretnej transformacji Fouriera                                #
###############################################################################

function dtft(freq::Real; signal::AbstractVector, fs::Real)
    dtft_val::ComplexF64 = 0.0
    for n in eachindex(signal)
        dtft_val += signal[n] * cispi(-2 * freq * n / fs)
    end
    return dtft_val
end


function dft(x::Vector)::Vector
    N = length(x)
    X = 1/N .* [sum(x[n + 1] * exp(-im * 2 * π * k * n / N) for n in 0:N-1) for k in 0:N-1]
    return X
end

function idft(x::Vector)::Vector
    N = length(x)
    y = 1/N .* [sum(x[k + 1] * exp(im * 2 * π * k * n / N) for k in 0:N-1) for n in 0:N-1]
    return y
end

function goertzel(x::Vector, k::Integer)::Complex
    N = length(x)
    ω = 2π * k / N
    coeff = 2 * cos(ω)
    s_p = 0.0
    s_pp = 0.0

    for n in 1:N
        s = x[n] + coeff * s_p - s_pp
        s_pp = s_p
        s_p = s
    end

    X = s_p - exp(-im * ω) * s_pp
    return X
end

function recursive_dft(N)
    signal = Complex{Float64}[]
    function f(x_n)
        push!(signal, x_n)
        current_N = length(signal)

        if current_N == 1
            return [sum(signal)]
        else
            X = Complex{Float64}(undef, current_N)
            for k in 1:current_N
                X[k] = sum([signal[n] * exp(-2im * π * (n-1) * (k-1) / current_N) for n in 1:current_N])
            end
            return X
        end
    end
    return f
end

function exp_recursive_dft(N, α)
    signal = Complex{Float64}[]
    w = exp(1im * 2 * π * α)

    function f(x_n)
        push!(signal, x_n)
        current_N = length(signal)

        if current_N == 1
            return [sum(signal)]
        else
            X = Complex{Float64}(undef, current_N)
            for k in 1:current_N
                X[k] = sum([signal[n] * w^(n-1) * exp(-2im * π * (n-1) * (k-1) / current_N) for n in 1:current_N])
            end
            return X
        end
    end
    return f
end

function cos_recursive_dft(N, α, a_m)
    signal = Complex{Float64}[]
    w = exp(1im * 2 * π * α)

    function f(x_n)
        push!(signal, x_n)
        current_N = length(signal)

        if current_N == 1
            return [sum(signal)]
        else
            X = Complex{Float64}(undef, current_N)
            for k in 1:current_N
                X[k] = sum([signal[n] * w^(n-1) * exp(-2im * π * (n-1) * (k-1) / current_N) for n in 1:current_N])
            end
            return X
        end
    end
    return f
end

function rdft(x::Vector)::Vector
    N = length(x)
    ζ = OffsetArray(
        [cispi(-2 * n / N) for n in 0:(N-1)],
        0:(N-1)
    )
    [
        sum((
            x[n+1] * ζ[(n*f)%N]
            for n in 0:(N-1)
        ))
        for f in 0:(N÷2)
    ]
end

function irdft(X::Vector, N::Integer)::Vector
    S = length(X)
    X1 = [n <= S ? X[n] : conj(X[2S-n+(N % 2 == 0 ? 0 : 1)]) for n in 1:N]
    real.(idft(X1))
end

function fft_radix2_dit_r(x::Vector)::Vector
        N = length(x)
    if N == 0 || (N & (N - 1)) != 0
        throw(ArgumentError("Length of input must be a power of 2"))
    end

    bits = Int(log2(N))
    function bitreverse(n, bits)
        reversed = 0
        for i in 1:bits
            reversed <<= 1
            reversed |= (n & 1)
            n >>= 1
        end
        return reversed
    end

    for i in 1:N
        j = bitreverse(i - 1, bits) + 1
        if i < j
            x[i], x[j] = x[j], x[i]
        end
    end

    m = 2
    while m <= N
        half_m = div(m, 2)
        w_m = cispi(-2 / m)
        for k in 1:m:N
            w = one(Complex{T})
            for j in 0:half_m-1
                t = w * x[k+j+half_m]
                u = x[k+j]
                x[k+j] = u + t
                x[k+j+half_m] = u - t
                w *= w_m
            end
        end
        m *= 2
    end

    return x
end

function ifft_radix2_dif_r(x::Vector)::Vector
    N = length(X)
    if N == 0 || (N & (N - 1)) != 0
        throw(ArgumentError("Length of input must be a power of 2"))
    end

    m = N
    while m >= 2
        half_m = div(m, 2)
        w_m = cispi(2 / m)
        for k in 1:m:N
            w = one(Complex{T})
            for j in 0:half_m-1
                t = X[k+j]
                u = X[k+j+half_m]
                X[k+j] = t + u
                X[k+j+half_m] = (t - u) * w
                w *= w_m
            end
        end
        m = half_m
    end

    bits = Int(log2(N))
    function bitreverse(n, bits)
        reversed = 0
        for i in 1:bits
            reversed <<= 1
            reversed |= (n & 1)
            n >>= 1
        end
        return reversed
    end

    for i in 1:N
        j = bitreverse(i - 1, bits) + 1
        if i < j
            X[i], X[j] = X[j], X[i]
        end
    end

    scale = inv(T(N))
    for i in 1:N
        X[i] *= scale
    end

    return X
end
end

function fft(x::Vector)::Vector
    FFTW.fft(x)
end

function ifft(x::Vector)::Vector
    FFTW.ifft(x)
end

function rfft(x::Vector)::Vector
    FFTW.rfft(x)
end

function irfft(x::Vector, N::Integer)::Vector
    FFTW.irfft(x)
end


###############################################################################
# Analiza częstotliwościowa sygnałów dyskretnych                              #
###############################################################################

fftfreq(N::Integer, fs::Real)::Vector = [n * N / fs for n in 0:(N-1)]
rfftfreq(N::Integer, fs::Real)::Vector = [n * N / fs for n in 0:(div(N,2))]

amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = abs.(fft(x .* w)) / (length(x) * mean(w))
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = amplitude_spectrum(x, w) .^ 2
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = abs2.(fft(x .* w)) / (sum(abs2, w) * fs)


function welch(x::Vector{T}, w::Vector{T}=rect(length(x)), L::Integer=1, fs::Real=1.0) where T
    N = length(x)
    K = length(w)
    Sxx = zeros(Complex{T}, N)

    for n in 1:L
        segment = x[(n-1)K+1:nK]
        segment = segment .* w
        Sxx += abs(fft(segment)) .^ 2
    end
    return Sxx / L
end

# Modelowanie sygnałów niestacjonarnych
function chirp_lin(t, f0, f1, T)
    return sin.(2 * π * (f0 * t .+ (f1 - f0) * t.^2 / (2 * T)))
end

function chirp_exp(t, f0, f1, T, φ)
    return sin.(2 * π * f0 * t .+ φ) .* exp.((f1 - f0) * t / T)
end

function stft(x::Vector{T}, w::Vector{T}, L::Integer) where T
    N = length(x)
    K = length(w)
    Sxx = Complex{T}(undef, K, L)

    for n in 1:L
        segment = x[(n-1)K+1:nK]
        segment = segment .* w
        Sxx[:, n] = rfft(segment)
    end

    return Sxx
end

function istft(X::AbstractMatrix{Complex{T}}, w::AbstractVector{T}, L::Integer=0, N::Integer=100)

    W_len = length(w)
    T_len = size(X, 2)
    y = zeros(Complex{T}, T_len * L + W_len)

    for i in 1:N
        for t in 1:T_len
            segment = X[:, t] .* exp(1im * angle(y[t * L - (W_len - 1) : t * L]))
            y[t * L - (W_len - 1) : t * L] .= y[t * L - (W_len - 1) : t * L] .+ ifft(segment)
        end
    end

    overlap_factor = sum(w.^2)
    y /= overlap_factor

    return real(y)
end


"""
Inputs:
    * X - spectrogram
    * w - STFT window
    * L - STFT overlap
    * N - number of Griffin-Lim iterations

Outputs:
    * Y - estimated STFT representation (amplitude + phase)
    * loss_values - a vector with reconsturction loss values (optional, recommended for dev)
"""
function phase_retrieval(X, w, L, N)
    ##
end


###############################################################################
# Systemy dyskretne                                                           #
###############################################################################

function lti_amp(f::Real, b::Vector, a::Vector)::Real
    M = length(b)
    K = length(a)
    num = sum(b[m+1] * cispi(-2 * f * m) for m in 0:M-1)
    denom = sum(a[k+1] * cispi(-2 * f * k) for k in 0:K-1)
    H_f = num / denom
    return abs(H_f)
end

function lti_phase(f::Real, b::Vector, a::Vector)::Real
    M = length(b)
    K = length(a)
    num = sum(b[m+1] * cispi(-2 * f * m) for m in 0:M-1)
    denom = sum(a[k+1] * cispi(-2 * f * k) for k in 0:K-1)
    H_f = num / denom
    return angle(H_f)
end

function conv(f::AbstractVector, g::AbstractVector)::Vector
    n = length(f)
    m = length(g)
    y = zeros(eltype(f), n + m - 1)
    for i in 1:n
        for j in 1:m
            y[i+j-1] += f[i] * g[j]
        end
    end
    return y
end

function fast_conv(f::Vector, g::Vector)::Vector
    N = length(f) + length(g) - 1

    f_padded = vcat(f, zeros(N - length(f)))
    g_padded = vcat(g, zeros(N - length(g)))

    F = fft(f_padded)
    G = fft(g_padded)
    Y = F .* G
    y = real(ifft(Y))

    return y
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    M = length(h)
    N = L + M - 1

    padded_h = vcat(h, zeros(N - M))
    H = fft(padded_h)

    y = zeros(eltype(x), length(x) + M - 1)

    for k in 1:L:length(x)
        xk = x[k:min(k + L - 1, end)]
        padded_xk = vcat(xk, zeros(N - length(xk)))
        Xk = fft(padded_xk)
        Yk = ifft(H .* Xk)
        yk_start = k
        yk_end = min(k + N - 1, length(y))
        y[yk_start:yk_end] += real(Yk[1:(yk_end-yk_start+1)])
    end

    return y[1:length(x)+M-1]
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    M = length(h)
    N = L + M - 1

    padded_h = vcat(h, zeros(N - M))
    H = fft(padded_h)

    y = []
    padded_x = vcat(zeros(M - 1), x, zeros(N - 1))

    for k in 1:L:(length(padded_x)-N+1)
        xk = padded_x[k:k+N-1]
        Xk = fft(xk)
        Yk = ifft(H .* Xk)
        y = vcat(y, real(Yk[M:end]))
    end

    return y[1:(length(x)+M-1)]
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    N = length(x)
    M = length(b) - 1
    K = length(a) - 1
    y = zeros(Float64, N)

    for n in 1:N
        for k in 0:M
            if n - k > 0
                y[n] += b[k+1] * x[n-k]
            end
        end
        for k in 1:K
            if n - k > 0
                y[n] -= a[k+1] * y[n-k]
            end
        end
    end
    return y
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    y_fwd_filt = lti_filter(b, a, x)
    y_rvs = reverse(y_fwd_filt)
    y_rvs_filt = lti_filter(b, a, y_rvs)
    return reverse(y_rvs_filt)
end

###############################################################################
# Projektowanie filtrów                                                       #
###############################################################################

function firwin_lp_I(order::Integer, F0::Float64)::Vector
    return [2F0 * sinc(2F0 * n) for n in -order/2:order/2]
end

function firwin_hp_I(order::Integer, F0::Float64)::Vector
    return [kronecker(Int(n)) - 2F0 * sinc(2F0 * n) for n in -order/2:order/2]
end

function firwin_bp_I(order::Integer, F1::Float64, F2::Float64)::Vector
    return [2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n) for n in -order/2:order/2]
end

function firwin_bs_I(order::Integer, F1::Float64, F2::Float64)::Vector
    return [kronecker(Int(n)) - (2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n)) for n in -order/2:order/2]
end

function firwin_lp_II(order::Integer, F0::Float64)::Vector
    N = range(start=order / 2, stop=order / 2, length=order)
    return [2F0 * sinc(2F0 * n) for n in N]
end

function firwin_bp_II(order::Integer, F1::Float64, F2::Float64)::Vector
    N = range(start=order / 2, stop=order / 2, length=order)
    return [2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n) for n in N]
end

firwin_diff(order::Int) = [n == 0 ? 0 : cospi(n) ./ n for n in -order:1:order]

function fir_I_LS(N::Integer, freq::Vector, amp::Vector, w::Vector)::Vector
    ##
end

function fir_I_remez(M::Integer, A::Function, W::Function)
    ##
end

function firwin_lp_kaiser_lp(f0, Δf, δ)
    ##
end