#= Zadanie 3:
#* correct solution
Dany jest idealny równomierny 5-bitowy kwantyzator q(x), którego najmniejszy poziom kwantyzacji ma
wartość a = -0.88, natomiast największy poziom kwantyzacji ma wartość b = 1.1. Oblicz sygnał błęd
kwantyzacji tego przetwornika dla dyskretnego sygnału x ∈ R^78. Jako rozwiazanie podaj energię sygnału
błędu kwantyzacji
=#

using CairoMakie

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


    kronecker(n) = ifelse(n == 0,1,0)
    triang_window(M) = [1 - abs(n)/(M + 1) for n in -M:M] 
        
    function fir_lp(
    order::Int = 56,
    fp::Float64 = 195.0,
    f0::Float64 = 25.35,
    z::Vector{Int} = [45, 23, 29, 15, 24, 15],
    )
        fn = f0/fp
        M = Int(order/2)
        h = [kronecker(n) - 2*fn*sinc(2π*fn*n) for n in -M:M]
        w = triang_window(M)
        hw = h .* w
        return hw
    end

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
    hw = fir_lp()
    HW = myDFT(hw)
    lines(abs.(HW))