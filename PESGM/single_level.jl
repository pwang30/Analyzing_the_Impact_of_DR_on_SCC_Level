# PESGM 2026, Canada 
# 19.Jul.2025 

import Pkg
using JuMP, BilevelJuMP,Plots,MAT, Ipopt

#-----------------------------------Define Parameters for Optimization-----------------------------------  
      
Load_total=[18.42,17.95,18.29,18.51,18.13,17.88,19.46,21.97,23.17,23.87,
23.91,23.77,23.80,23.82,24.23,23.79,26.01,26.91,25.26,23.69,22.12,20.04,18.17,18.01]*10   # (MW)
T=length(Load_total)

Pˢᴳₘₐₓ=[6.584, 5.760, 3.781, 3.335, 3.252, 2.880]*15            
Pˢᴳₘᵢₙ=[3.292, 2.880, 1.512, 0.667, 0.650, 0.288]*15                        
Oᵐ=[6.20, 7.10, 10.47, 12.28, 13.53, 15.36]                        


#-----------------------------------Define Model-----------------------------------
model= Model()

#-------Define Primal Variales

@variable(model, 0<=Pˢᴳ¹[1:T]<=Pˢᴳₘₐₓ[1])    
@variable(model, 0<=Pˢᴳ²[1:T]<=Pˢᴳₘₐₓ[2])    
@variable(model, 0<=Pˢᴳ³[1:T]<=Pˢᴳₘₐₓ[3])                   
@variable(model, 0<=Pˢᴳ⁴[1:T]<=Pˢᴳₘₐₓ[4])                
@variable(model, 0<=Pˢᴳ⁵[1:T]<=Pˢᴳₘₐₓ[5])            
@variable(model, 0<=Pˢᴳ⁶[1:T]<=Pˢᴳₘₐₓ[6])
               

#-------Define Primal Constraints

Power_balance=Dict()
for t in 1:T
    Power_balance[t]=@constraint(model, Pˢᴳ¹[t]+ Pˢᴳ²[t] +Pˢᴳ³[t] +Pˢᴳ⁴[t] +Pˢᴳ⁵[t] +Pˢᴳ⁶[t]== Load_total[t])     # power balance , dual variable: λᴱₜ                             
@constraint(model, Pˢᴳ¹[t].<=Pˢᴳₘₐₓ[1])           # bounds for the output of SGs 
@constraint(model, Pˢᴳₘᵢₙ[1].<=Pˢᴳ¹[t]) 
@constraint(model, Pˢᴳ²[t].<=Pˢᴳₘₐₓ[2])           
@constraint(model, Pˢᴳₘᵢₙ[2].<=Pˢᴳ²[t])           
@constraint(model, Pˢᴳ³[t].<=Pˢᴳₘₐₓ[3])       
@constraint(model, Pˢᴳₘᵢₙ[3].<=Pˢᴳ³[t])         
@constraint(model, Pˢᴳ⁴[t].<=Pˢᴳₘₐₓ[4])       
@constraint(model, Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁴[t])         
@constraint(model, Pˢᴳ⁵[t].<=Pˢᴳₘₐₓ[5])       
@constraint(model, Pˢᴳₘᵢₙ[5].<=Pˢᴳ⁵[t])        
@constraint(model, Pˢᴳ⁶[t].<=Pˢᴳₘₐₓ[6])       
@constraint(model, Pˢᴳₘᵢₙ[6].<=Pˢᴳ⁶[t])
end

#-------Define Objective Functions 
@variable(model, operation_cost[1:T])

for t in 1:T
@constraint(model, operation_cost[t]== Oᵐ[1]*Pˢᴳ¹[t] +Oᵐ[2]*Pˢᴳ²[t] +Oᵐ[3]*Pˢᴳ³[t] +Oᵐ[4]*Pˢᴳ⁴[t] +Oᵐ[5]*Pˢᴳ⁵[t] +Oᵐ[6]*Pˢᴳ⁶[t])
end

obj_Primal=sum(operation_cost)

@objective(model, Min, obj_Primal)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Ipopt.Optimizer)
optimize!(model)


price_original=zeros(1,T)
for t in 1:T
    price_original[t]=dual(Power_balance[t])
end

plot(price_original')
plot(Load_total)



matwrite("Load_total.mat", Dict("Load_total" => Load_total))


