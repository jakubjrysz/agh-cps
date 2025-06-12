using Random

function sawtooth_wave(t)
    return -2*((t+0.5)-floor(t+0.5)-0.5)
end

t = 0:0.1:3

println(sawtooth_wave.(t))

