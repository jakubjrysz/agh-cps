using CairoMakie
using LinearAlgebra


function rozwiazanie(;
    b::Vector{Float64} = [0.0424461855453743, -0.2546771132722458, 0.6366927831806145, -0.848923710907486, 0.6366927831806145, -0.2546771132722458, 0.0424461855453743],
    a::Vector{Float64} = [1.0, -0.5376591778088213, 1.3966595537904305, 0.012955823584018755, 0.7085675481324825, 0.09935532488214605, 0.3066620341017809],
    F::Vector{Float64} = [0.03, 0.17, 0.32, 0.44],
)
                                                    #z = exp(j2π*f)     gdzie f = f/fp 
    f_len = length(F)
    b_len = length(b)
    a_len = length(a)
    licznik = zeros(ComplexF64, f_len)
    mianownik = ones(ComplexF64, f_len)
    gain = zeros(f_len)

    for i in 1:f_len
        z = exp(im * 2 * π * F[i])
        for j in 1:b_len
            licznik[i] += b[j] * z ^ (- j + 1)
        end
        for k in 2:a_len
            mianownik[i] += a[k] * z ^ (- k + 1)
        end
        gain[i] = abs(licznik[i]/mianownik[i])
    end

    avg_gain = sum(gain) / f_len

    return avg_gain
end


rozwiazanie()