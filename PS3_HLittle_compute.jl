#I've modeled this compute file off of those we've seen in previous problem states

using Parameters, Plots
include("PS3_HLittle_model.jl")

prim, res = Initialize_R()
#distrib = Initialize_M()
price = Initialize_P()
#kgrid = Initialize_K()

##### Test ######
println(prim.J)

#=
#originally thought I should use the BI_Solve function but I don't need to do that ...?
for i in prim.J:-1:1
    a=i
    Bellman(prim, res, price; a) #very likely have a problem here...semicolon ??
end #close for loop
=#

Bellman(prim, res, price) #very likely have a problem here...semicolon ??

#just unpacking what I need for now
@unpack val_func = res
@unpack kap = kgrid

Plots.plot(kap, val_func[:, :, 50])