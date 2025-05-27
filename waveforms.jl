using CairoMakie
using LinearAlgebra


ramp_wave(t::Real)::Real = 2 * (t % 1)

triangle_wave(t::Real)::Real = 4*(0.5 - abs((t % 1) - 0.5))-1

sawtooth_wave(t::Real)::Real = 2 * (t % 1) - 1

