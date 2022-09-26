#I've modeled this compute file off of those we've seen in previous problem states

using Parameters, Plots
include("PS3_HLittle_model.jl")

prim, res = Initialize_R()
#distrib = Initialize_M()
price = Initialize_P()
#kgrid = Initialize_K()

##### Test ######
println(prim.J)

Bellman(prim, res, price) #keep having problems here

#just unpacking what I need for now
@unpack val_func = res
@unpack kap = prim

Plots.plot(kap, val_func[:, :, 65])