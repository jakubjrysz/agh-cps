using CairoMakie

function fftfreq(N,fp)
    frequencyResolution = fp/N
    frequencyVector =[0:frequencyResolution:fp]

    return frequencyVector
end

X = fftfreq(2000,10)