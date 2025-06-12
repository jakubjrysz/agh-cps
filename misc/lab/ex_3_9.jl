using Random

function ramp_wave(t)
    if t%1==0
        return 0 
    else
        return t - floor(t)
    end
end

t = 0:0.1:3

println(ramp_wave.(t))

