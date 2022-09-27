#I've modeled this compute file off of those we've seen in previous problem states

using Parameters, Plots
include("PS3_HLittle_model.jl")

prim, res = Initialize_R()
price = Initialize_P()
distrib = Initialize_M()

Bellman(prim, res, price) #I was indexing with k prime when I should have been using k

Fill_Mu(prim, res, distrib) #function to fill out the cross sectional distribution


#unpacking what I'll use to plot
@unpack val_func, pol_func  = res
@unpack kap, nk = prim
@unpack mass = distrib

#Plot the value function at age 50
Plots.plot(kap, val_func[:, :, 50], labels=["High State" "Low State"])

#Plot the savings rate at age 20
disc_k_grid = zeros(nk, 2)
disc_k_grid[:, 1] = kap
disc_k_grid[:, 2] = kap
disc_k_grid = 0.94.*disc_k_grid
savings20 = pol_func[:,:,20] - disc_k_grid
Plots.plot(kap, savings20, labels=["High State" "Low State"])

#Plot the policy function at age 20
Plots.plot(kap, pol_func[:,:,20], labels=["High State" "Low State"])




