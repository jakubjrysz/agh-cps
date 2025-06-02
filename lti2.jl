using LinearAlgebra

function r_roznicowe(;
    b::Vector{Float64} = [0.08221128890062884, -0.25714826784473266, 0.3586135586310597, -0.2571482678447327, 0.08221128890062886],
    a::Vector{Float64} = [1.0, 0.4830157588163749, 0.5907268075519141, 0.09596213126178924, 0.025583754648032666],
    x::Vector{Float64} = [0.8, -0.76, 0.55, -0.06, -0.35, -0.03, 0.31, -0.05, 0.61, -0.71, -0.09, 0.26, -0.05, -0.88, -0.69, 0.58, -0.64, 0.23, -0.38, -0.41, -0.96, -0.5, -0.21, 0.49, 0.53, -0.62, -0.35, -0.88, -0.79, 0.24, 0.15, -0.63, -0.67, -0.66, 0.18, -0.54, -0.03, 0.26, 0.69, -0.36, 0.57, -0.43, -0.75, -0.01, 0.54, -0.69, -0.62, -0.81, 0.29],
    L::Int = 66,
)
    M = length(b)
    K = length(a)
    N = length(x)
    y = zeros(N)

    for n in 1:N
        for m in 1:M
            if n - m + 1 > 0
                y[n] += b[m] * x[n - m + 1]
            end
        end
        
        for k in 2:K
            if n - k + 1 > 0 
                y[n] -= a[k] * y[n - k + 1]
            end
        end
    end

    return y
end

function zera_bieguny(;
    zz::Vector{ComplexF64} = ComplexF64[-1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.4486778207816606 + 0.8471249785044375im, 0.4486778207816606 - 0.8471249785044375im, 0.5892575170333751 + 0.6488708837470728im, 0.5892575170333751 - 0.6488708837470728im, 0.7685166023732889 + 0.2574225101383183im, 0.7685166023732889 - 0.2574225101383183im],
    k::Float64 = 0.0008961796507509385,
    F::Vector{Float64} = [0.01, 0.03, 0.32],
)
z_len = length(zz)
p_len = length(pp)
f_len = length(F)
licz = ones(ComplexF64, f_len)
mian = ones(ComplexF64, f_len)

for i in 1:f_len
    z = exp(-im*2π*F[i])
    licz[i] = zz[1]
    for k in 2:z_len
        licz[i] *=  (1 - zz[k]*z)
    end
    for k in 1:p_len
        mian[i] *= (1 - pp[k]*z)
    end
end

return licz[1]

end

zera_bieguny()