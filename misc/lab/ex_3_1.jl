using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

fs = 1000
dt = 1/fs
t = 0.25:dt:0.25+255*dt
x = 2*sin.(50*pi*t .+ pi/4)

println(length(x))
println(x)

scatter(t, x)