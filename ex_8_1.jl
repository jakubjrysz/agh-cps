using CairoMakie

#=
frequency resolution w DFT jest dane wzorem:

f = 1/(deltaT*N)

N - liczba punktów DFT
deltaT - odstęp między próbkami w czasie
=#

delta_t = 50e-06
N = 8192

f = 1/(delta_t*N)   #Hz