using CairoMakie
using LinearAlgebra


    ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
    sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
    triangular_wave(t::Real)::Real = ifelse(mod(t + 1 / 4, 1.0) < 1 / 2, 4mod(t + 1 / 4, 1.0) - 1, -4mod(t + 1 / 4, 1.0) + 3)
    square_wave(t::Real)::Real = ifelse(mod(t, 1) < 0.5, 1, -1)

    t = 0:0.01:2

    scatter(t,sawtooth_wave.(t))