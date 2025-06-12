using CairoMakie

delta_t = 0.001
sampled_length = 2000

f_res = 1/(delta_t*sampled_length)

println(f_res)

a = f_res * 0 
b = f_res * 200
c = f_res * 1000
d = f_res * 1600