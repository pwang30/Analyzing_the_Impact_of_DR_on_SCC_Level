# PESGM 2026, Canada 
# 19.Jul.2025 

import Pkg
#Xusing JuMP,Gurobi, Bilevel, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT, CPLEX, Ipopt
using JuMP, BilevelJuMP,Plots,MAT, Gurobi


#----------------------------------Parameters----------------------------------#
load_original=[18.42,17.95,18.29,18.51,18.13,17.88,19.46,21.97,23.17,23.87,
23.91,23.77,23.80,23.82,24.23,23.79,26.01,26.91,25.26,23.69,22.12,20.04,18.17,18.01]*10   # (MW)

T=length(Load_total)

Pˢᴳₘₐₓ=[6.584, 5.760, 3.781, 3.335, 3.252, 2.880]*15            
Pˢᴳₘᵢₙ=[3.292, 2.880, 1.512, 0.667, 0.650, 0.288]*15                        
Oᵐ=[6.20, 7.10, 10.47, 12.28, 13.53, 15.36]                         

b_max=1.2
b_min=0.8
# Define TOU period type for each hour: 1=valley, 2=flat, 3=peak
tou_type = [
    1,1,1,1,1,1,                # 1–6h: valley
    2,2,2,                      # 7–9h: flat
    3,3,3,3,3,3,3,3,3,3,3,      # 10–20h: peak
    2,2,                        # 21–22h: flat
    1,1                         # 23–24h: valley
]


self_elasticity = -0.25       # own-price elasticity
cross_elasticity_adj = 0.05   # adjacent-period elasticity
cross_elasticity_nonadj = 0.02 # non-adjacent weak elasticity

E = zeros(Float64, T, T) # Construct Elasticity Matrix (24x24)
for i in 1:T
    for j in 1:T
        if i == j
            E[i, j] = self_elasticity
        elseif abs(i - j) == 1
            E[i, j] = cross_elasticity_adj
        else
            E[i, j] = cross_elasticity_nonadj
        end
    end
end

#-----------------------------------Define Bilevel Model-----------------------------------
model= BilevelModel()


#---Lower level---
@variable(Lower(model), Pˢᴳ¹[1:T])    
@variable(Lower(model), Pˢᴳ²[1:T])    
@variable(Lower(model), Pˢᴳ³[1:T])                   
@variable(Lower(model), Pˢᴳ⁴[1:T])                
@variable(Lower(model), Pˢᴳ⁵[1:T])            
@variable(Lower(model), Pˢᴳ⁶[1:T])
            
Power_balance=Dict()
for t in 1:T
    Power_balance[t]=@constraint(Lower(model), Pˢᴳ¹[t]+ Pˢᴳ²[t] +Pˢᴳ³[t] +Pˢᴳ⁴[t] +Pˢᴳ⁵[t] +Pˢᴳ⁶[t]== P_D[t])     # power balance , dual variable: λᴱₜ                             
@constraint(Lower(model), Pˢᴳ¹[t].<=Pˢᴳₘₐₓ[1])           # bounds for the output of SGs 
@constraint(Lower(model), Pˢᴳₘᵢₙ[1].<=Pˢᴳ¹[t]) 
@constraint(Lower(model), Pˢᴳ²[t].<=Pˢᴳₘₐₓ[2])           
@constraint(Lower(model), Pˢᴳₘᵢₙ[2].<=Pˢᴳ²[t])           
@constraint(Lower(model), Pˢᴳ³[t].<=Pˢᴳₘₐₓ[3])       
@constraint(Lower(model), Pˢᴳₘᵢₙ[3].<=Pˢᴳ³[t])         
@constraint(Lower(model), Pˢᴳ⁴[t].<=Pˢᴳₘₐₓ[4])       
@constraint(Lower(model), Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁴[t])         
@constraint(Lower(model), Pˢᴳ⁵[t].<=Pˢᴳₘₐₓ[5])       
@constraint(Lower(model), Pˢᴳₘᵢₙ[5].<=Pˢᴳ⁵[t])        
@constraint(Lower(model), Pˢᴳ⁶[t].<=Pˢᴳₘₐₓ[6])       
@constraint(Lower(model), Pˢᴳₘᵢₙ[6].<=Pˢᴳ⁶[t])
end

@objective(Lower(model), Min, sum(Oᵐ[1].*Pˢᴳ¹) +sum(Oᵐ[2].*Pˢᴳ²) +sum(Oᵐ[3].*Pˢᴳ³) +sum(Oᵐ[4].*Pˢᴳ⁴) +sum(Oᵐ[5].*Pˢᴳ⁵) +sum(Oᵐ[6].*Pˢᴳ⁶)  )  # obj of lower-level problem: min generation cost


#---Upper level---

new_load = base_load .* (1 .+ ΔL)


@variable(Upper(model), Δp[1:T])
@variable(Upper(model), q[1:T])
@variable(Upper(model), UL_obj[1:T])

for t in 1:T
    @constraint(Upper(model), UL_obj[t]==λ[t]*P_tou[t])
    @constraint(Upper(model), Δp[t]== (λ[t]-price_original[t])/price_original[t] )   # new price - original price
    @constraint(Upper(model), q[t]== load_original[t] *(1 + E[t,:]*Δp ) )               # new load - original load
    @constraint(Upper(model), b_min*q[t] <=  P_D[t] <= b_max*q[t]  )                       # new load

end

@objective(Upper(model), Min, sum(UL_obj)  )                # obj of lower-level problem: min energy cost


#---------------------------solve the bi-level model-------------------

BilevelJuMP.set_mode(model, BilevelJuMP.StrongDualityMode())
set_optimizer(model, Gurobi.Optimizer)
optimize!(model)


#---------------------------Get the results-------------------
Pˢᴳ²=JuMP.value.(Pˢᴳ²)
Pˢᴳ³=JuMP.value.(Pˢᴳ³)
Pˢᴳ⁴=JuMP.value.(Pˢᴳ⁴)
Pˢᴳ⁵=JuMP.value.(Pˢᴳ⁵)
Pˢᴳ²⁷=JuMP.value.(Pˢᴳ²⁷)
Pˢᴳ³⁰=JuMP.value.(Pˢᴳ³⁰)
k=JuMP.value.(k)

P_D=JuMP.value.(P_D)


plot!(P_D)

matwrite("SCL_prices_bus30_disp_MAC.mat", Dict("SCL_prices_bus30_disp_MAC" => SCL_prices_bus30_MAC))
