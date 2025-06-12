using Makie
using LinearAlgebra
using GLMakie
using CairoMakie

fs = 2048
dt = 1/fs

t0 = 5
t1 = 10

t = t0:dt:t1
x = 0.25*sin.(pi/2*t .+ pi)

length(x)