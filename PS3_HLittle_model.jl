##Note that this problem set builds closely off of PS4 from 712 
##I refer to the matlab code on the Canvas website for 712

using Parameters, Plots #import the libraries we want

#Using this space to think on μ
μ_holder = zeros(66, 1)
μ_holder[1]=1
for i = 2:66
    μ_holder[i]=μ_holder[i-1]/(1.011)
end #end for loop
μ_sum = sum(μ_holder)
μ_holder = μ_holder*(1/μ_sum)
println(sum(μ_holder)) #sanity check

#Create a primitives struct
@with_kw struct Primitives
    #demographics
    J::Int64 = 66 #lifespan
    JR::Int64 = 46 #age of retirement
    tR::Int64 = J-JR+1 #length of retirement
    tW::Int64 = JR-1 #total 
    age_vec::Array{Int64, 1} = collect(1:1:66)
    ret_vec::Array{Int64, 1} = collect(46:1:65)
    work_vec::Array{Int64, 1} = collect(1:1:45)
    n::Float64 = 0.011 #population growth rate

    #preferences
    β::Float64 = 0.97 #discount factor
    σ::Float64 = 2 #coeff of rel risk aversion
    γ::Float64 = 0.42 #weight on consumption, note this is gamma not y

    #tax
    θ::Float64 = 0.11 #tax to fund social security benefits

    #for capital law of motion
    α::Float64 = 0.36 #capital share, see cobb douglas production function
    δ::Float64 = 0.06 #rate of depreciation

    #for high and low states
    hs::Float64 = 3 #high state productivity
    ls::Float64 = 0.5 #low state producitivity

    #state transitions, markov transition probabilities decomposed
    Phh::Float64 = 0.9261
    Phl::Float64 = 0.0739
    Pll::Float64 = 0.9811
    Plh::Float64 = 0.0189

    #age efficiency profile
    e::Array{Float64, 1} = [      0.59923239 
    0.63885106 
    0.67846973 
    0.71808840 
    0.75699959 
    0.79591079 
    0.83482198 
    0.87373318 
    0.91264437 
    0.95155556 
    0.99046676 
    0.99872065 
     1.0069745 
     1.0152284 
     1.0234823 
     1.0317362 
     1.0399901
     1.0482440 
     1.0564979 
     1.0647518 
     1.0730057 
     1.0787834 
     1.0845611 
     1.0903388 
     1.0961165 
     1.1018943 
     1.1076720 
     1.1134497 
     1.1192274 
     1.1250052 
     1.1307829 
     1.1233544 
     1.1159259 
     1.1084974 
     1.1010689 
     1.0936404 
     1.0862119 
     1.0787834 
     1.0713549 
     1.0639264
     1.0519200
     1.0430000
     1.0363000
     1.0200000
     1.0110000]

    #capital grid
    maxkap::Float64 = 14 #max value in capital grid, should not be binding! be careful across experiments
    minkap::Float64 = 0 #minimum value of capital grid, changed matlab code from 0.01 to 0 because of what problem set asks
    nk::Int64 = 180 #number of grid points
    kap::Array{Float64, 1} = collect(minkap:((maxkap-minkap)/(nk-1)):maxkap) #this has a uniform density of points, unlike matlab code

    #initialize a very small number
    neg::Float64 = -1e10

    #Keeping track of the weights
    μ_vec = μ_holder

end #end primitives struct

#Create a mutable struct for results
mutable struct Results
    val_func::Array{Float64, 3} #value function, note that this "3" means we're initializing a 3 dimensional array 
    pol_func::Array{Float64, 3} #policy function, 3 dimensions
    lab_func::Array{Float64, 3} #labor supply function, 3 dimensions
end #close mutable struct

#Create a mutable struct for the distribution
mutable struct Distribution
    mass::Array{Float64, 3} #cross sectional distribution, also 3 dimensions
end #close mutable struct

#Create a mutable struct for the prices
mutable struct Prices
    K_agg::Float64 #aggregate capital
    L_agg::Float64 #aggregate labor
    w::Float64 #wages
    r::Float64 #interest rate
    b::Float64 #benefit amount
end #close mutable struct


function Initialize_R()
prim = Primitives()
val_func = zeros(prim.nk, 2, prim.J) #three dimensions of zeros, capital grid by states by number of ages
pol_func = zeros(prim.nk, 2, prim.J) #like above
lab_func = zeros(prim.nk, 2, prim.J) #like above
res = Results(val_func, pol_func, lab_func)
return prim, res #note that we're returning both primitives and res here so running Initialize_R also calls our primitives
end #close function

function Initialize_P()
    prim = Primitives() 
    K_agg = 3.639 #initial guess for benchmark SS
    L_agg = 0.4308 #initial guess for benchmark SS
    #below, we use the aggregate levels of capital and labor 
    w = (1-prim.α)*K_agg^prim.α*L_agg^(-prim.α)
    rk = prim.α*K_agg^(prim.α-1)*L_agg^(1-prim.α)
    r = rk - prim.δ #net off depreciation, rk is paid by firms, r is seen by hh
    ret_mass = sum(prim.μ_vec[46:66])
    b = (prim.θ*w*L_agg)/(ret_mass)
    price = Prices(K_agg, L_agg, w, r, b)
    return price
end #close function

function Initialize_M()
    prim = Primitives() 
    #initialize the whole thing as zeros
    mass = zeros(prim.nk, 2, prim.J)
    distrib = Distribution(mass) #initialize results struct, did have q at one point
    return distrib #return distribution
end #close function

#Garrett's Get Index function
function get_index(val::Float64, grid::Array{Float64,1})
    n = length(grid)
    index = 0 #preallocation
    if val<=grid[1] #LEQ smallest element 
        index = 1
    elseif val>=grid[n] #GEQ biggest element
        index = n
    else
        index_upper = findfirst(z->z>val, grid)
        index_lower = index_upper - 1
        val_upper, val_lower = grid[index_upper], grid[index_lower] #values
        index = index_lower + (val - val_lower)  / (val_upper - val_lower) #weighted average
    end
    index #return
end #close the get_index function

function Fill_Mu(prim::Primitives,res::Results, distrib::Distribution) #note that this must accout for the generations growing over time...

    #initialize the first generation
    distrib.mass[1, 1, 1] = 0.2037 #about twenty percent begin with zero assets and high prod
    distrib.mass[1, 2, 1] = 0.7963 #about 8 percent begin with zero assets and low prod
    for p_index = 2:66
        #focus on those coming from the high state
        pol_func_high = res.pol_func[:, 1, p_index-1] #create a vector of the capital choices in high state the previous age
        for pos_h = 1:prim.nk #looking over all choices of capital by INDEX
            kprime = pol_func_high[pos_h] #gives the actual choice of k prime at each index from yesterday
            pos = round(Int64, get_index(kprime, prim.kap)) #maps the k prime choice from yesterday into an index of the capital grid
            println(pos)
            #the high prod people at a given capital level yesterday (pos_h) go to high state today (pos) at k prime choice with prob Phh
            distrib.mass[pos, 1, p_index] += prim.Phh*distrib.mass[pos_h, 1, p_index-1] 
            #the people at a given capital level today via HL come from the mass of those who chose it from high state last period
            distrib.mass[pos, 2, p_index] += prim.Phl*distrib.mass[pos_h, 1, p_index-1] 
        end #end loop over the capital holdings

        #focus on those coming from the low state
        pol_func_low = res.pol_func[:, 2, p_index-1] #create a vector of the capital choices in low state at specific age
        for pos_l = 1:prim.nk
            kprime = pol_func_low[pos_l] #which value of k prime are you looking for
            pos = round(Int64, get_index(kprime, prim.kap))             
            #logic starts on the RHS, bring the low to high people to where their policy func puts them (LHS)
            distrib.mass[pos, 1, p_index] += prim.Plh*distrib.mass[pos_l, 2, p_index-1] 
            #logic starts on the RHS, bring the low to low people to where their policy func puts them (LHS)
            distrib.mass[pos, 2, p_index] += prim.Pll*distrib.mass[pos_l, 2, p_index-1] 
        end #end loop over the capital holdings
    
    end #end loop over age panels

    
    #now reweight each of the age panels
    for i = 1:66
        weight = prim.μ_vec[i]
        distrib.mass[:, :, i] = distrib.mass[:, :, i].*weight
    end #end rewight loop
    

    return distrib.mass
end #close function

#############################################################################################
#find the optimal policy functions and value functions
#the following section follows the backward iteration we saw in lecture 8 of the computational bootcamp

#Bellman Operator
function Bellman(prim::Primitives,res::Results, price::Prices) #note that a is going to be an age profile, we'll iterate over it later...NOT accounting for age here
    @unpack val_func, pol_func, lab_func = res #unpack value function
    @unpack w, r, b = price #unpack the interest rate
    @unpack kap, nk, J, JR, tR, tW, β, σ, γ, θ, α, δ, hs, ls, Phh, Phl, Pll, Plh, e, age_vec= prim #unpack model primitives

        ###
        #Solve the problem of the final period people
        ###
            #we know that k'=0 for every single one of these people
            for i = 1:nk
                #k_hold = kap[i] #kept getting an error when I tried to put this into c directly below
                c = b + (1+r)*kap[i] #consumption
                val_func[i, :, J] .= c^((1-σ)γ)/(1-σ) #the value function (for both the high and low productivity states) is simply the utility given the capital you're coming in with
                #val_func[i, 2, J] = c^((1-σ)*γ)/(1-σ)
            end #end for loop

        ###
        #Solve the problem of the retired people, age 46 to 65
        ###
        for a = 65:-1:46 #loop over retired but not 66
            for k = 1:nk #for every capital value today
                candidate_max = -Inf #bad candidate max (states do not matter here)
                budget = b + (1+r)*kap[k]
                for kp = 1:nk
                    #alt approach, should we have depreciation in the household's BC/capital law of motion
                    c = budget - kap[kp] #consumption associated with a given choice of k prime
                    if c > 0 #only consider this if consumption if greater than zero
                        val = (c^((1-σ)γ))/(1-σ) + β*val_func[kp, 1, a+1] #flow utility plus value func associated with kp choice at tomorrow's age (not productivity state dependent here, so I pick first col)
                        if val > candidate_max
                            val_func[k, :, a] .= val #record value associated with that choice
                            pol_func[k, :, a] .= kap[kp] #record the choice
                            candidate_max = val #update candidate max
                        end #close innner if statement
                    end #close if statement
                end #close for loop, looping over k prime
            end #close for loop, looping over k
        end #close for loop over retired but not 66
        ###
        #Solve the problem of the working age people
        ###
        for a = 45:-1:1 #loop over working age
            for k = 1:nk #for every capital value today
                candidate_max_h = -Inf #bad candidate max from high state
                candidate_max_l = -Inf #bad candidate max from low state
                
                #focusing only on the high productivity shock
                for kp_h = 1:nk
                    l = (γ*(1-θ)hs*e[a]*w - (1-γ)*((1+r)*kap[k]-kap[kp_h]))/((1-θ)*w*hs*e[a])
                    if l>1
                        l = 1
                    elseif l<0
                        l = 0
                    end #end if loop
                    budget_h = (1-θ)l*hs*e[a]*w + (1+r)*kap[k]
                    #alt approach, should we have depreciation in the houshold's BC/capital law of motion
                    c = budget_h - kap[kp_h] #consumption associated with a given choice of k prime (post tax lab inc, cap inc, less cap invest)
                    if c > 0 #only consider this if consumption if greater than zero
                        val = ((c^γ*(1-l)^(1-γ))^(1-σ))/(1-σ) + β*Phh*val_func[kp_h, 1, a+1] + β*Phl*val_func[kp_h, 2, a+1] #flow util plus discounted expectation of val tomorrow given k prime choice
                        if val > candidate_max_h
                            val_func[k, 1, a] = val #record value associated with that choice, just the high state now
                            pol_func[k, 1, a] = kap[kp_h] #record the choice, just the high state now
                            lab_func[k, 1, a] = l #less sure on this, if all else is optimal here, then that's what you're choosing for labor too
                            candidate_max_h = val #update candidate max
                        end #close innner if statement
                    end #close if statement
                end #close for loop, looping over k prime
                
                #focusing only on the low producitivity shock
                for kp_l = 1:nk
                    l = (γ*(1-θ)*ls*e[a]*w - (1-γ)*((1+r)*kap[k]-kap[kp_l]))/((1-θ)*w*ls*e[a])
                    if l>1
                        l = 1
                    elseif l<0
                        l = 0
                    end #end if loop
                    #alt approach, should we have depreciation in the houshold's BC/capital law of motion
                    budget_l = (1-θ)l*ls*e[a]*w + (1+r)*kap[k]
                    c = budget_l - kap[kp_l] #consumption associated with a given choice of k prime (post tax lab inc, cap inc, less cap invest)
                    if c > 0 #only consider this if consumption if greater than zero
                        val = ((c^γ*(1-l)^(1-γ))^(1-σ))/(1-σ) + β*Plh*val_func[kp_l, 1, a+1] + β*Pll*val_func[kp_l, 2, a+1] #flow util plus discounted expectation of val tomorrow given k prime choice
                        if val > candidate_max_l
                            val_func[k, 2, a] = val #record value associated with that choice, just the high state now
                            pol_func[k, 2, a] = kap[kp_l] #record the choice, just the high state now
                            lab_func[k, 2, a] = l #less sure on this, if all else is optimal here, then that's what you're choosing for labor too
                            candidate_max_l = val #update candidate max
                        end #close innner if statement
                    end #close if statement
                end #close for loop, looping over k prime
            end #close for loop, looping over k
        end #close for loop over working age
end #end the Bellman Operator


##############################################################################





