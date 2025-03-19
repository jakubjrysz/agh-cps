using Random
using CairoMakie
using LinearAlgebra


function SNR(p_signal,p_noise)
    return 10*log10(p_signal/p_noise)
end

SNR(1000,1)