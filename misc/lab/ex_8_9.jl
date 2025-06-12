using CairoMakie

function fftfreq(N,fp)
    frequencyResolution = fp/N
    frequencyVector =[0:frequencyResolution:fp/2]

    return frequencyVector
end

X = fftfreq(2000,10)