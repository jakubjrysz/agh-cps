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
    #hw = fir_lp()
    #HW = myDFT(hw)
    #lines(abs.(HW))
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
function poly_from_roots(r::AbstractVector)
        p = [1.0 + 0im]
        for i in eachindex(r)
            p = conv(p, [1, -r[i]])
        end
        return p
    end

        function rozwiazanie_7_2(;
        zz::Vector{ComplexF64}=ComplexF64[0.9562125977134925-0.2926729710342487im, 0.6431481557756844+0.7657417643842708im, 0.9562125977134925+0.2926729710342487im, 0.6431481557756844-0.7657417643842708im, 0.9494540691139096-0.3139059901356441im, 0.6842289643189376+0.7292672516896902im, 0.9494540691139096+0.3139059901356441im, 0.6842289643189376-0.7292672516896902im, 0.9154123917349581-0.4025172705090848im, 0.8016754572775892+0.5977595345946131im, 0.9154123917349581+0.4025172705090848im, 0.8016754572775892-0.5977595345946131im],
        pp::Vector{ComplexF64}=ComplexF64[0.6064983645211166+0.787230370569453im, 0.9585153738231527-0.277256566212053im, 0.6064983645211166-0.787230370569453im, 0.9585153738231527+0.277256566212053im, 0.5363569616798629+0.7946801478013233im, 0.9537353148226411-0.25407497984321314im, 0.5363569616798629-0.7946801478013233im, 0.9537353148226411+0.25407497984321314im, 0.12166088338468982+0.6657516938860896im, 0.9290620344841986-0.14298335418965522im, 0.12166088338468982-0.6657516938860896im, 0.9290620344841986+0.14298335418965522im],
        k::Float64=0.28119079368088185,
        F::Vector{Float64}=[0.14, 0.17, 0.2],
    )
        a = poly(pp)
        b = poly(zz) .* k
        return abs(sum(a)/sum(b))
        end

        rozwiazanie_7_2()