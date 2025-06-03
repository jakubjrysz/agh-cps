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
#=
ZAD1.

Dany jest dyskretny sygnał x € C39, którego próbki zostały probrane z ciągłego sygnału o 
ograniczonym paśmie, z szybkością fp = 858 próbek na sekunde.
 Oblicz 39-punktową dyskretną transformację Fouriera tego sygnału oraz znajdź wartości dyskretnego widma 
Fouriera tego sygnału dla następujących częstotliwość f = [-242, 44, 396, 418]. Jako rozwiązanie podaj sumę faz tych
składowych częstotliwościowych.

DZIAŁA
=#

function zad1(;
    fp::Int = 693,
    x::Vector{ComplexF64} = ComplexF64[-0.65 + 0.13im, -1.02 + 0.82im, -0.9 + 0.88im, -1.49 + 0.52im, 0.68 + 0.37im, -0.16 - 0.52im, 0.37 + 1.18im, 0.0 - 0.77im, 0.31 + 0.09im, -0.25 - 0.38im, -0.6 - 0.21im, -0.66 - 0.23im, -0.25 - 0.1im, 0.56 + 0.25im, -0.07 + 0.04im, -1.35 - 0.18im, 0.02 - 0.8im, 0.56 + 0.22im, 0.82 + 1.47im, 0.63 + 1.05im, -0.61 + 0.47im, -0.8 - 0.17im, -0.5 - 0.87im, 0.88 + 0.2im, -0.35 + 0.09im, -0.79 + 0.02im, 0.48 + 0.51im, 0.17 - 0.48im, 0.2 - 1.1im, 0.26 - 0.9im, 0.19 + 1.27im, -0.0 + 0.78im, -1.2 - 0.36im],
    f::Vector{Int} = [-189, -168, -42, 105, 252, 273],
)
    X = myDFT(x)
    N = length(x)
    k = [dft_freq_to_index(f[i], fp, N) for i in eachindex(f)]
    phase = 0
    for i in k
        phase += angle(X[i])
    end
    return phase
end
#=
ZAD2.

Oblicz średnią dyskretnego sygnału z R370. 
Dyskretny sygnał z powstał w wyniki pobrania N = 471 próbek z ciągłego sygnału 
y(t) = 2.0 · g(2.9 · t - 4.3) z szybkością fp = 457.14 próbek na sekundę.
Pierwsza próbka 1 = y(t1) została pobrana w chwili t1 = 7.9.
Funkcja g(t) zwraca wartości sygnału bipolarnej fali prostokątnej o następujących parametrach:
amplituda 1, okres 1 sekunda, składowa stała 0, 9(0) = 0, oraz g(t) = 1 dla t € (0, 1).

DZIAŁA
=#

function zad2(;
    fp::Float64 = 457.14,
    t1::Float64 = 7.9,
    N::Int = 471,
)
    function g(t)
        t_mod = mod(t , 1)
        return ifelse(t_mod == 0, 0, t_mod < 0.5 ? 1 : -1) 
    end

    t = t1 .+ (0:N-1) ./ fp
    y = 2.0 .* g.(2.9 .* t .- 4.3)

    return mean(y)
end


#=
ZAD3.

Dany jest dyskretny system liniowy niezmienny w czasie, który jest 
opisany Z-transmitacją H(z). 
Transmitancja H(z) jest zdefiniowana poprzez dwa wektory b = [bo, b1,..., bм] € RM+1 oraz
a = [a0, a1,...,N] ERN+1,
które są odpowiednio współczynnikami wielomianu 
w liczniku i mianowniku funkcji wymiernej H(z). Oblicz przesunięcie fazowe (f) tego systemu 
dla częstotliwości ƒER znormalizowanej względem częstotliwości 
próbkowania. Jako odpowiedź podaj średnie przesunięcie fazowe 
dla następujących częstotliwości F = [0.02, 0.36, 0.38], to znaczy,

DZIAŁA
=#

function zad3(;
    b::Vector{Float64} = [0.024733551701640596, 0.011781887468929471, 0.05249319238023622, 0.027475783356817555, 0.05249319238023622, 0.011781887468929471, 0.024733551701640603],
    a::Vector{Float64} = [1.0, -3.0220628341265607, 5.209220125491564, -5.523324466370808, 3.867298543047256, -1.6562250598356836, 0.35566068215132135],
    F::Vector{Float64} = [0.02, 0.36, 0.38],
)
    faza = lti_phase(b,a,F)
    return mean(faza)
end


#=
ZAD4.

Dany jest idealny równomierny 4-bitowy kwantyzator q(x),
którego najmniejszy poziom kwantyzacji ma wartość a = 0.0077,
natomist największy poziom kwantyzacji ma wartość b = 1.0. 
Oblicz sygnał błęd kwantyzacji tego przetwornika dla dyskretnego sygnału x R54. 
Jako rozwiązanie podaj wartość skuteczną sygnału błędu kwantyzacji.

DZIAŁA (pamiętać o zmianie N - ilosci bitów)
=#

function zad4(;
    a::Float64 = -2.7,
    b::Float64 = 3.8,
    x::Vector{Float64} = [2.876, 3.08473, 3.29347, 3.5022, 3.71094, -2.68033, -2.47159, -2.26286, -2.05412, -1.84539, -1.63665, -1.42792, -1.21918, -1.01045, -0.80171, -0.59298, -0.38424, -0.17551, 0.03323, 0.24196, 0.4507, 0.65943, 0.86817, 1.0769, 1.28564, 1.49437, 1.70311, 1.91184, 2.12058, 2.32931, 2.53805, 2.74678, 2.95552, 3.16425, 3.37299, 3.58172, 3.79046, -2.60081, -2.39207, -2.18334, -1.9746, -1.76587, -1.55713, -1.3484, -1.13966, -0.93093, -0.72219, -0.51346, -0.30472, -0.09599, 0.11275, 0.32148, 0.53022, 0.73895, 0.94769],
)
    error = quantization_error(a, b, 4, x)
    return energy(error)
end

#=

ZAD5.

Oblicz odpowiedź impulsową hЄ R7 nierekursywnego filtru pasmowoprzepustowego rzędu 76 
o liniowej charakterystyce fazowej. Filtr zaprojektuj tak aby przy częstotliwości próbkowania fp = 197.0 Hz, 
3 dB pasmo przepustowe było między częstotliwościami f1 = 3.94 i f2 = 51.22 Hz. Do zaprojektowania filtru wykorzystaj metodę okien czasowych i okno Hanninga. 
Jako rozwiązanie podaj sumę wartości wektora ho indeksach z = [27, 21, 20, 10].

=#

function zad5(;
    order::Int = 76,
    fp::Float64 = 197.0,
    f1::Float64 = 3.94,
    f2::Float64 = 51.22,
    z::Vector{Int} = [27, 21, 20, 10],
)
    M = div(order, 2)
    h = fir_bp(fp, f1, f2, order)
    w = hann_window(M)
    hw = h .* w
    sum = 0
    for i in z
        sum += hw[i]
    end
    return sum
end

zad5()