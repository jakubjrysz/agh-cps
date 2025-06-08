using CairoMakie
using LinearAlgebra


ramp_wave(t::Real)::Real = 2 * (t % 1)

triangle_wave(t::Real)::Real = abs((t % 1) - 0.5) 

sawtooth_wave(t::Real)::Real = 2 * (t % 1) - 1

t = 0:0.01:2

scatter(t,triangle_wave.(t))