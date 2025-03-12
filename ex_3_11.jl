using GLMakie

function triangular_wave(t)
    return 2*abs(2*((4*t+0.25) - floor(0.75+4*t)))-1
end

x = 0:0.0001:1

println(triangular_wave.(x))

lines(triangular_wave.(x))