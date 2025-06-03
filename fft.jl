function zad5(;
    order::Int = 77,
    fp::Float64 = 197.0,
    f1::Float64 = 3.94,
    f2::Float64 = 51.22,
    z::Vector{Int} = [27, 21, 20, 10],
)
    M = div(order, 2)  # środek filtru
    h = fir_bp(fp, f1, f2, order)
    w = hann_window(M)
    hw = h .* w

    # Zmieniamy indeksy z [-M:M] na [1:order] (Julia indeksuje od 1)
    sum_hw = sum(hw[i] for i in z)
    return sum_hw
end

# Okno Hanninga o symetrii względem zera: n ∈ -M:M
hann_window(M) = [0.5 - 0.5*cos(2π * n / (2M + 1)) for n in -M:M]

# Filtr pasmowoprzepustowy: odejmujemy dwa filtry dolnoprzepustowe
function fir_bp(fp, f1, f2, order)
    fn1 = f1 / fp
    fn2 = f2 / fp
    M = div(order, 2)
    return [2fn2 * sinc(2fn2 * n) - 2fn1 * sinc(2fn1 * n) for n in -M:M]
end

zad5()