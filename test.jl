#= Zadanie 3:
#* correct solution
Dany jest idealny równomierny 5-bitowy kwantyzator q(x), którego najmniejszy poziom kwantyzacji ma
wartość a = -0.88, natomiast największy poziom kwantyzacji ma wartość b = 1.1. Oblicz sygnał błęd
kwantyzacji tego przetwornika dla dyskretnego sygnału x ∈ R^78. Jako rozwiazanie podaj energię sygnału
błędu kwantyzacji
=#


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