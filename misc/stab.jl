using LinearAlgebra
using CairoMakie

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
stab()