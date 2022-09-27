##Note that this problem set builds closely off of PS4 from 712 
##I refer to the matlab code on the Canvas website for 712

using Parameters, Plots #import the libraries we want

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
val_func = zeros(180, 2, 66) #three dimensions of zeros, capital grid by states by number of ages
pol_func = zeros(180, 2, 66) #like above
lab_func = zeros(180, 2, 66) #like above
res = Results(val_func, pol_func, lab_func)
return prim, res #note that we're returning both primitives and res here so running Initialize_R also calls our primitives
end #close function

function Initialize_P()
    K_agg =
    L_agg = 
    
    w = 1.05
    r = 0.05
    b = 0.2
    price = Prices(K_agg, L_agg, w, r, b)
    return price
end #close function

function Initialize_M() #note that this must accout for the generations growing over time...
    prim = Primitives() #initialize primtiives
    for i = 1:prim.J
        #first, put a uniform distribution (over states and capital) for all ages
        mass[:, :, i] = ones(prim.nk, 2) 
        mass[:, :, i] = mass[:, :, i]./(2*prim.nk)
        #then, adjust for ages using the given rate of population growth
        ###Actually, don't think I need to do this here...
    end #close for loop
    distrib = Distribution(mass) #initialize results struct, did have q at one point
    distrib #return deliverables
end

#############################################################################################
#find the optimal policy functions and value functions
#the following section follows the backward iteration we saw in lecture 8 of the computational bootcamp

#Bellman Operator
function Bellman(prim::Primitives,res::Results, price::Prices) #note that a is going to be an age profile, we'll iterate over it later...NOT accounting for age here
    @unpack val_func, pol_func, lab_func = res #unpack value function
    @unpack w, r, b = price #unpack the interest rate
    @unpack kap, nk, J, JR, tR, tW, β, σ, γ, θ, α, δ, hs, ls, Phh, Phl, Pll, Plh, e, age_vec= prim #unpack model primitives
    #v_next = zeros(nk, 2, J) #next guess of value function to fill, accounting for each state of the world

    #for a in reverse(age_vec) #big loop across all ages

        #if a == J #if you're in the terminal age
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

        #elseif  J > a >= JR #if you're retired
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
        #else a < JR
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
        #end #end else if ... got rid of this approach
    #end #close the loop over ages ... got rid of this approach
end #end the Bellman Operator


##############################################################################





