using GLMakie

function square_wave(t)
    return 2*(2*floor(t) - floor(2*t))+1
end

x = 0:0.0001:1

println(square_wave.(x))

lines(square_wave.(x))
