using Makie
using LinearAlgebra

mean(x::AbstractVector)::Real = sum(x)/length(x)

energy(x::AbstractVector)::Real = sum(abs2, x)

power(x::AbstractVector)::Real = energy(x)/length(x)

rms(x::AbstractVector)::Real = sqrt(power(x))




ramp_wave(t::Real)::Real = 2 * (t % 1)

triangle_wave(t::Real)::Real = 4*(0.5 - abs((mod(t,1)) - 0.5)) - 1

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
    amp = zeros(f_len)

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

function lti_phase(b::Vector, a::Vector, F::Vector)
                                                    #z = exp(j2π*f)     gdzie f = f/fp 
    f_len = length(F)
    b_len = length(b)
    a_len = length(a)
    licznik = zeros(ComplexF64, f_len)
    mianownik = ones(ComplexF64, f_len)
    phase = zeros(f_len)

    for i in 1:f_len
        z = exp(im * 2 * π * F[i])
        for j in 1:b_len
            licznik[i] += b[j] * z ^ (- j + 1)
        end
        for k in 2:a_len
            mianownik[i] += a[k] * z ^ (- k + 1)
        end
        phase[i] = angle(licznik[i]/mianownik[i])
    end
    return phase
end

function convolution(x, h)
    n = length(x)
    m = length(h)
    y = zeros(n + m - 1)
    for i in 1:n
        for j in 1:m
            y[i + j - 1] += x[i] * h[j]
        end
    end
    return y
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

function quantization_error(a, b, N, x)
    L = range(start = a, stop = b, length = 2 ^ N)
    x_quantized = [quantize(xi, L) for xi in x]
    error = x .- x_quantized
    return error
end

function dft_freq(N, fs)
    k = collect(-div(N, 2) : div(N-1, 2))
    return fs * k / N
end

function dft_freq_to_index(freq, fs, N)
    k = round(Int, freq * N / fs)
    return mod(k, N) + 1
end



kronecker(M) = ifelse(n == 0,1,0)
tri_window(M) = [1 - abs(n)/(M+1) for n in -M:M]
hann_window(M) = [0.5 - 0.5*cos(2*π*n/(2*M+1)) for n in -M:M]
hamm_window(M) = [0.54 - 0.46*cos(2*π*n/(2*M+1)) for n in -M:M]
blackman_window(M) = [0.42 + 0.5*cos(2*π*n/(2*M+1)) + 0.08*cos(4*π*n/(2*M+1)) for n in -M:M]


function fir_lp(fp, f0, order)
    fn = f0/fp
    M = div(order, 2)
    h = [2*fn * sinc(2*π*fn*n) for n in -M:M]
    return h
end

function fir_hp(fp, f0, order)
    fn = f0/fp
    M = div(order, 2)
    h = [kronecker(n) - 2*fn * sinc(2*π*fn*n) for n in -M:M]
    return h
end

function fir_bp(fp,f1,f2,order)
    fn1 = f1/fp
    fn2 = f2/fp
    M = div(order, 2)
    h = [2*fn2*sinc(2*fn2*n) - 2*fn1*sinc(2*fn1*n) for n in -M:M]
    return h
end

function fir_bs(fp,f1,f2,order)
    fn1 = f1/fp
    fn2 = f2/fp
    M = div(order, 2)
    h = [kronecker(n) - (2*fn2*sinc(2*π*fn2*n) - 2*fn1*sinc(2*π*fn1*n)) for n in -M:M]
    return h
end

function rozwiazanie32(;
    b::Vector{Float64} = [0.6668547023844128, -1.720764471909391, 3.4806618241912046, -3.8658925116288825, 3.4806618241912037, -1.720764471909391, 0.6668547023844127],
    a::Vector{Float64} = [1.0, -2.2575274288819127, 3.948876910415752, -3.813445803637163, 2.951472093696976, -1.2364482229285918, 0.3946840490385092],
    F::Vector{Float64} = [0.14, 0.43, 0.48, 0.5],
)
    b_len = length(b)
    a_len = length(a)
    f_len = length(F)
    mianowniki = zeros(ComplexF64, f_len)
    liczniki = zeros(ComplexF64, f_len)
    y = zeros(ComplexF64, f_len)

    for i in 1:f_len
        z = exp(im * 2 * π * F[i])
        for k in 1:b_len
            liczniki[i] += b[k] * z^(-k + 1)
        end
        for k in 1:b_len
            mianowniki[i] += a[k] * z^(-k + 1)
        end
        y[i] = liczniki[i]/mianowniki[i]
    end
    return sum(abs.(y))/f_len
end
#0.9475142159299259 dobrze

function rozwiazanie33(;
    zz::Vector{ComplexF64} = ComplexF64[0.5653138377631579 + 0.8248759087483948im, 0.5653138377631579 - 0.8248759087483948im, 0.9090424796534677 + 0.41670345593176317im, 0.9090424796534677 - 0.41670345593176317im],
    pp::Vector{ComplexF64} = ComplexF64[-0.27610745027671657 + 0.6519220973669482im, -0.27610745027671657 - 0.6519220973669482im, -0.1706544392027319 + 0.2080237764058478im, -0.1706544392027319 - 0.2080237764058478im],
    k::Float64 = 0.05804559733761661,
    F::Vector{Float64} = [0.01, 0.44, 0.46],
)
    z_len = length(zz)
    p_len = length(pp)
    f_len = length(F)
    mianowniki = ones(ComplexF64, f_len)
    liczniki = ones(ComplexF64, f_len)
    y = ones(ComplexF64, f_len)

    for i in 1:f_len
        z = exp(im * 2 * π * F[i])
        for k in 1:p_len
            mianowniki[i] *= (1 - pp[k] * z^(-1))
        end
        for k in 1:z_len
            liczniki[i] *= (1 - zz[k] * z^(-1))
        end
        y[i] = k*liczniki[i]/mianowniki[i]
    end
    return sum(abs.(y))/f_len
end
#0.6676898085110858 dobrze

function rozwiazanie34(;
    a::Float64 = -2.4,
    b::Float64 = 4.9,
    x::Vector{Float64} = [-1.84, -1.3681, -0.89621, -0.42431, 0.04758, 0.51948, 0.99137, 1.46327, 1.93516, 2.40706, 2.87895, 3.35085, 3.82275, 4.29464, 4.76654, 4.76157, 4.28967, 3.81778, 3.34588, 2.87399, 2.40209, 1.9302, 1.4583, 0.98641, 0.51451, 0.04261, -0.42928, -0.90118, -1.37307, -1.84497, -2.31686, -2.41124, -1.93935, -1.46745, -0.99556, -0.52366, -0.05176, 0.42013, 0.89203, 1.36392, 1.83582, 2.30771, 2.77961, 3.2515, 3.7234, 4.19529, 4.66719, 4.86092, 4.38902, 3.91712, 3.44523, 2.97333],
)
    L = range(start = a, stop = b, length = 2^5)
    function quantize(x, L)
        return L[argmin(abs.(L .- x))]
    end
    xq = [quantize(x[i],L) for i in eachindex(x)]
    err = x .- xq
    pwr = sum(err.^2)/length(err)
    return sqrt(pwr)
end
#0.07288866395547805 dobrze

function rozwiazanie35(;
    order::Int = 86,
    fp::Float64 = 192.0,
    f0::Float64 = 32.64,
    z::Vector{Int} = [26, 24, 82],
)
    M = div(order,2)
    fn = f0/fp
    h = [2 * fn * sinc(2 * fn * n) for n in -M:M]
    w = [0.42 + 0.5*cos((2 * π * n)/(2 * M + 1)) + 0.08*cos((4 * π * n)/(2 * M + 1)) for n in -M:M]
    hw = h .* w
    out = sum(hw[i] for i in z)
    return out
end
#0.006983941709695499 dobrze

function rozwiazanie36(;
    b::Vector{Float64} = [0.15220776723622292, -0.6306612471456423, 1.1561685097078895, -1.1561685097078895, 0.6306612471456425, -0.15220776723622292],
    a::Vector{Float64} = [1.0, -1.1397725742429585, 1.51371415095299, -0.3171939570028648, 0.16876931161974204, 0.2613749456390448],
    x::Vector{Float64} = [-0.12, -0.57, -0.18, -0.44, -0.71, -0.46, -0.38, 0.9, 0.43, 0.28, -0.41, 0.73, -0.7, -0.24, -0.63, -0.06, -0.85, 0.07, 0.76, -0.37, 0.99, 0.21, 0.26, -0.16, -0.01, -0.66, 0.48],
    L::Int = 54,
)
    b_len = length(b)
    a_len = length(a)
    N = length(x)
    y = zeros(L)

    for i in 1:L
        for j in 1:b_len
            if i - j + 1 > 0 && i - j + 1 <= N
                y[i] += b[j] * x[i - j + 1]
            end
        end
        for j in 2:a_len
            if i - j + 1 > 0
                y[i] -= a[j] * y[i - j + 1]
            end
        end
    end
    return sum(y.^2)
end
#3.3525141321526473 dobrze

function rozwiazanie37(;
    fp::Int = 893,
    x::Vector{ComplexF64} = ComplexF64[-0.05 + 0.15im, -0.51 - 0.06im, 0.46 + 0.01im, 0.9 + 0.31im, -0.56 + 0.38im, 0.87 + 0.61im, 0.53 + 0.4im, -1.63 + 1.23im, 0.31 + 0.77im, 1.43 + 0.06im, 0.43 + 0.66im, -1.15 + 0.41im, -0.47 - 0.39im, -1.86 - 0.44im, 1.38 - 0.15im, 0.18 + 0.08im, 0.1 - 1.8im, -0.9 - 0.57im, -0.04 + 1.66im, -0.97 + 1.17im, -0.51 - 1.03im, 0.15 - 0.8im, -0.09 - 0.5im, -0.9 + 1.25im, -0.2 + 1.12im, -0.24 + 1.15im, 1.38 + 0.03im, 0.53 - 0.76im, 0.46 + 0.58im, 0.59 + 0.89im, -1.09 + 0.51im, 0.08 + 0.89im, -0.39 - 0.49im, 0.53 + 0.71im, -0.59 - 0.26im, -0.42 - 0.3im, 0.14 - 0.52im, 0.85 - 0.05im, -0.21 - 0.58im, -0.1 + 0.32im, 0.13 - 1.0im, 1.26 + 0.05im, 0.01 + 0.26im, -0.2 - 1.52im, -1.42 + 0.38im, -0.21 - 1.1im, 0.35 - 0.54im],
    f::Vector{Int} = [133, 171, 323, 418],
)
    N = length(x)
    f_len = length(f)
    fn = [f[i]/fp for i in 1:f_len]

    X_all = [sum(x[n + 1] * exp(-im * 2 * π * k * n / N) for n in 0:N-1) for k in 0:N-1]
    X_spec = [sum(x[n + 1] * exp(-im * 2 * π * fi * n) for n in 0:N-1) for fi in fn]

    return sum(angle.(X_spec))

end
#-4.361081197373714 dobrze

function rozwiazanie38(;
    m::Vector{Float64} = [2.8, 2.8019, 2.8038, 2.8057, 2.8076, 2.8095, 2.8114, 2.8133, 2.8152, 2.8171, 2.819, 2.8209, 2.8228, 2.8247, 2.8266, 2.8285, 2.8304, 2.8323, 2.8342, 2.8361, 2.838, 2.8399, 2.8418, 2.8437, 2.8456, 2.8475, 2.8494, 2.8513, 2.8532, 2.8551, 2.857, 2.8589, 2.8608, 2.8627, 2.8646, 2.8665, 2.8684, 2.8703, 2.8722, 2.8741, 2.876, 2.8779, 2.8798, 2.8817, 2.8836, 2.8855, 2.8874, 2.8893, 2.8912, 2.8931, 2.895, 2.8969, 2.8988, 2.9007, 2.9026, 2.9045, 2.9064, 2.9083, 2.9102, 2.9121, 2.914],
    s::Vector{Float64} = [0.0025, 0.0588, 0.8666, 0.7365, 0.1205, 0.3544, 0.7307, 0.4551, 0.0726, 0.4234, 0.1429, 0.4319, 0.8926, 0.447, 0.9485, 0.1193, 0.5426, 0.824, 0.7288, 0.5842, 0.351, 0.5389, 0.0885, 0.5073, 0.4411, 0.2934, 0.5071, 0.5116, 0.0003, 0.9083, 0.7886, 0.3045, 0.8251, 0.8614, 0.617, 0.6531, 0.7969, 0.1681, 0.2543, 0.6885, 0.927, 0.451, 0.7601, 0.9688, 0.8076, 0.2111, 0.6989, 0.7451, 0.489, 0.0804, 0.1332, 0.4731, 0.2665, 0.6928, 0.5649, 0.9544, 0.3258, 0.4117, 0.3996, 0.7947, 0.0048],
    t::Vector{Float64} = [2.8133, 2.86023, 2.88664, 2.89443, 2.84085, 2.84807],
)
    delta_t = m[2] - m[1]

    x = [sum(s[i] * sinc((ti - m[i])/delta_t) for i in eachindex(s)) for ti in t]
    return sum(x)
end
#2.010477765167538 dobrze

function rozwiazanie39(;
    x::Vector{Float64} = [-3.7, -3.08, -0.03, -4.52, 1.96, -2.09, -0.57, 4.51, -2.29, -3.16, -0.33, 0.09, 4.28, -1.48, -2.41, 0.25, -4.49, -4.92, -2.37, -3.47, 4.74, 1.43, -2.83, 1.64, -3.8, -4.85, 0.3, -4.99, -3.69, -1.57, -0.28, -2.36, -4.77, 0.97, -0.35, -0.26, 0.61, -3.76, -3.19, 4.72, -2.3, 0.0, 1.58, -0.13, 2.75, -4.13, -3.22, -3.15, 2.5, -0.71, 3.49, -2.69, -2.04, -2.29, -2.43, 3.01, 2.17, 2.21, -4.16, -1.48, -0.69, -3.52, -3.6],
    h::Vector{Float64} = [-0.64, 4.88, -0.42, -1.42, 0.06, -3.67, -3.69, -2.31, 4.06, -1.02, -2.1],
)
    N = length(x)
    K = length(h)
    y = zeros(N + K - 1)

    for i in 1:N
        for j in 1:K
            y[i + j - 1] += x[i] * h[j]
        end
    end
    return sum(y)/length(y)
end
#5.9212849315068485 dobrze

function rozwiazanie40(;
    fp::Float64 = 471.46,
    t1::Float64 = 0.28,
    N::Int = 922,
)
    g(t) = ifelse(mod(t,1) < 0.5, 1, -1)
    t = range(start = t1, length = N, step = 1/fp)
    y = 3.0 .* g.(2.5 .* t .- 1.0)
    return sum(y)/length(y)
end
#0.07158351409978309 dobrze

function rozwiazanie41(;
    z::Vector{ComplexF64} = ComplexF64[-1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    p::Vector{ComplexF64} = ComplexF64[0.26243533330430674 + 0.693265210762646im, 0.26243533330430674 - 0.693265210762646im, 1.5669479036478788 + 2.5582560988402854im, 0.19233484328103675 - 0.3140128556269943im, 0.17452793889436508 + 0.0im],
    k::Float64 = 0.019847795197400566,
)
    p_abs = abs.(p)

    if any(p_abs .> 1)
        return -1.0
    else
        return 1.0
    end
end
#-1.0 dobrze

function rozwiazanie42(;
    b::Vector{Float64} = [0.0036995705131426322, -0.006588701086142232, 0.005123154642098176, -2.0536741824576278e-19, -0.005123154642098175, 0.006588701086142232, -0.003699570513142633],
    a::Vector{Float64} = [1.0, -2.978210743260216, 5.773741495194624, -6.608511452910479, 5.494231273099657, -2.69585928272244, 0.8614467267882088],
    F::Vector{Float64} = [0.22, 0.23, 0.24, 0.31, 0.32, 0.36],
)
    a_len = length(a)
    b_len = length(b)
    f_len = length(F)
    mi = zeros(ComplexF64, f_len)
    li = zeros(ComplexF64, f_len)
    y = zeros(ComplexF64, f_len)

    for i in 1:f_len
        z = exp(im * 2 * π * F[i])
        for j in 1:b_len
            li[i] += b[j] * z^(-j + 1)
        end
        for j in 1:a_len
            mi[i] += a[j] * z^(-j + 1)
        end
        y[i] = li[i]/mi[i]
    end
    return sum(abs.(y))/f_len
end
#0.002486328200146685 dobrze

#zadanie ze stabilnością ab

function stab(;
        b::Vector{Float64}=[0.20496999142745434, -1.0248499571372718, 2.0496999142745436, -2.0496999142745436, 1.0248499571372718, -0.20496999142745434],
        a::Vector{Float64}=[1.0, -2.04328699669025, 2.205686901762555, -1.0925050400038385, 0.2795267139181644, 0.06196592669626961],
    )

    #a = mianownik, b = licznik
    N = length(a) - 1;
    companion_matrix = zeros(N,N)
    a_reversed = reverse(a)

    for i in 1:N-1
        companion_matrix[i+1, i] = 1
    end
    for i in 1:N
        companion_matrix[i, N] = -a_reversed[i]
    end

    #display(companion_matrix)
    p = eigen(companion_matrix)
    roots = abs.(p.values)

    if any(roots .> 1)
        return -1.0
    else
        return 1.0
    end
end