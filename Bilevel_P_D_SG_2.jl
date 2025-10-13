# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Now, this is the economic_dispatch+UC WITHOUT SCL constraints on a modified IEEE-30 bus system, SCL on each bus is calculated in an offline process according to the solution
# 20.March.2025

import Pkg
using JuMP,Gurobi,Plots,MAT



#-----------------------------------Define Parameters for Optimization-----------------------------------  
      
Load_total=[18.42,17.95,18.29,18.51,18.13,17.88,19.46,21.97,23.17,23.87,
23.91,23.77,23.80,23.82,24.23,23.79,26.01,26.91,25.26,23.69,22.12,20.04,18.17,18.01]*10^3/8   # (MW)
T=length(Load_total)

Pˢᴳₘₐₓ=[5.584, 4.760, 3.781, 2.335]*10^3/4             # Max generation of SGs                    
Pˢᴳₘᵢₙ=[3.292, 2.880, 1.512, 0.667]*10^3/5           #  Min generation of SGs                   
Kˢᵗ=[200, 125, 92.5, 72]*10^3 /10                           #  Startup cost of SGs                    
Kˢʰ=[50, 28.5, 18.5, 14.4]*10^3  /10                        #  Shutdown cost of SGs                   
Oᵐ=[6.20, 7.10, 10.47, 12.28]                        #  Marginal generation cost of SG   
Oⁿˡ=[17.431, 15.005, 13.755, 10.930]*10^3 /10         #  No-load cost of SGs    
yˢᴳ₀=[0, 0, 0, 0]                                       # Initial status of SGs (1-startup, 0-shutdown)               
bᵐᵃˣ= 2                                                            # Bidding variable_energy upper bound



#-----------------------------------Define Primal-Dual Model-----------------------------------
model = Model(Gurobi.Optimizer)

#-------Define Primal Variales
@variable(model, b[1:T])          

@variable(model, Pˢᴳ¹[1:T])        # generation of SGs 
@variable(model, Pˢᴳ²[1:T])  
@variable(model, Pˢᴳ³[1:T])                     
@variable(model, Pˢᴳ⁴[1:T])              

@variable(model, yˢᴳ¹[1:T],Bin)      # status of SGs
@variable(model, yˢᴳ²[1:T],Bin)          
@variable(model, yˢᴳ³[1:T],Bin)           
@variable(model, yˢᴳ⁴[1:T],Bin)                                  

@variable(model, Cᵁ¹[1:T]>=0)                 # startup costs and shutdown costs for SGs                
@variable(model, Cᴰ¹[1:T]>=0)                                 
@variable(model, Cᵁ²[1:T]>=0)                 
@variable(model, Cᴰ²[1:T]>=0)  
@variable(model, Cᵁ³[1:T]>=0)   
@variable(model, Cᴰ³[1:T]>=0)
@variable(model, Cᵁ⁴[1:T]>=0)
@variable(model, Cᴰ⁴[1:T]>=0)        
          


#-------Define Dual Variales
@variable(model, λᴱ[1:T]>=0)            # price for market clearing                

@variable(model, μᵐⁱⁿˢᴳ¹[1:T]>=0)
@variable(model, μᵐᵃˣˢᴳ¹[1:T]>=0)
@variable(model, μᵐⁱⁿˢᴳ²[1:T]>=0)
@variable(model, μᵐᵃˣˢᴳ²[1:T]>=0)
@variable(model, μᵐⁱⁿˢᴳ³[1:T]>=0)
@variable(model, μᵐᵃˣˢᴳ³[1:T]>=0)
@variable(model, μᵐⁱⁿˢᴳ⁴[1:T]>=0)
@variable(model, μᵐᵃˣˢᴳ⁴[1:T]>=0)

@variable(model, σˢᵗˢᴳ¹[1:T]>=0)
@variable(model, σˢʰˢᴳ¹[1:T]>=0)
@variable(model, σˢᵗˢᴳ²[1:T]>=0)
@variable(model, σˢʰˢᴳ²[1:T]>=0)
@variable(model, σˢᵗˢᴳ³[1:T]>=0)
@variable(model, σˢʰˢᴳ³[1:T]>=0)
@variable(model, σˢᵗˢᴳ⁴[1:T]>=0)
@variable(model, σˢʰˢᴳ⁴[1:T]>=0)

@variable(model, ψᵐᵃˣˢᴳ¹[1:T]>=0)
@variable(model, ψᵐᵃˣˢᴳ²[1:T]>=0)
@variable(model, ψᵐᵃˣˢᴳ³[1:T]>=0)
@variable(model, ψᵐᵃˣˢᴳ⁴[1:T]>=0)



#-------Define Primal Constraints
@constraint(model, b.>=1)                    # bound for bidding decicions
@constraint(model, b.<=bᵐᵃˣ)

@constraint(model, Pˢᴳ¹+Pˢᴳ²+Pˢᴳ³+Pˢᴳ⁴ == Load_total)     # power balance 

@constraint(model, Pˢᴳ¹.<=yˢᴳ¹*Pˢᴳₘₐₓ[1])           # bounds for the output of SGs with UC 
@constraint(model, yˢᴳ¹*Pˢᴳₘᵢₙ[1].<=Pˢᴳ¹)
@constraint(model, Pˢᴳ².<=yˢᴳ²*Pˢᴳₘₐₓ[2])       
@constraint(model, yˢᴳ²*Pˢᴳₘᵢₙ[2].<=Pˢᴳ²)             
@constraint(model, Pˢᴳ³.<=yˢᴳ³*Pˢᴳₘₐₓ[3])       
@constraint(model, yˢᴳ³*Pˢᴳₘᵢₙ[3].<=Pˢᴳ³)       
@constraint(model, Pˢᴳ⁴.<=yˢᴳ⁴*Pˢᴳₘₐₓ[4])       
@constraint(model, yˢᴳ⁴*Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁴)     

@constraint(model, Cᵁ¹[1]>=(yˢᴳ¹[1]-yˢᴳ₀[1])*Kˢᵗ[1])        # startup costs and shutdown costs for SGs , dual variables: σˢᵗₜ , σˢʰₜ
@constraint(model, Cᴰ¹[1]>=(yˢᴳ₀[1]-yˢᴳ¹[1])*Kˢʰ[1])  
@constraint(model, Cᵁ²[1]>=(yˢᴳ²[1]-yˢᴳ₀[1])*Kˢᵗ[2])        
@constraint(model, Cᴰ²[1]>=(yˢᴳ₀[1]-yˢᴳ²[1])*Kˢʰ[2]) 
@constraint(model, Cᵁ³[1]>=(yˢᴳ³[1]-yˢᴳ₀[1])*Kˢᵗ[3])
@constraint(model, Cᴰ³[1]>=(yˢᴳ₀[1]-yˢᴳ³[1])*Kˢʰ[3])
@constraint(model, Cᵁ⁴[1]>=(yˢᴳ⁴[1]-yˢᴳ₀[1])*Kˢᵗ[4])
@constraint(model, Cᴰ⁴[1]>=(yˢᴳ₀[1]-yˢᴳ⁴[1])*Kˢʰ[4])

for t in 2:T
    @constraint(model, Cᵁ¹[t]>=(yˢᴳ¹[t]-yˢᴳ¹[t-1])*Kˢᵗ[1])        
    @constraint(model, Cᴰ¹[t]>=(yˢᴳ¹[t-1]-yˢᴳ¹[t])*Kˢʰ[1])  
    @constraint(model, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*Kˢᵗ[2])        
    @constraint(model, Cᴰ²[t]>=(yˢᴳ²[t-1]-yˢᴳ²[t])*Kˢʰ[2]) 
    @constraint(model, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*Kˢᵗ[3])
    @constraint(model, Cᴰ³[t]>=(yˢᴳ³[t-1]-yˢᴳ³[t])*Kˢʰ[3])
    @constraint(model, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*Kˢᵗ[4])
    @constraint(model, Cᴰ⁴[t]>=(yˢᴳ⁴[t-1]-yˢᴳ⁴[t])*Kˢʰ[4])
end



#-------Define Dual Constraints   
@constraint(model, Oⁿˡ[1] -Pˢᴳₘₐₓ[1]*μᵐᵃˣˢᴳ¹[T] +Pˢᴳₘᵢₙ[1]*μᵐⁱⁿˢᴳ¹[T] +Kˢᵗ[1]*σˢᵗˢᴳ¹[T] -Kˢʰ[1]*σˢʰˢᴳ¹[T]+ ψᵐᵃˣˢᴳ¹[T] >=0)             # dual constraints for UC, when t==T
@constraint(model, Oⁿˡ[2] -Pˢᴳₘₐₓ[2]*μᵐᵃˣˢᴳ²[T] +Pˢᴳₘᵢₙ[2]*μᵐⁱⁿˢᴳ²[T] +Kˢᵗ[2]*σˢᵗˢᴳ²[T] -Kˢʰ[2]*σˢʰˢᴳ²[T]+ ψᵐᵃˣˢᴳ²[T] >=0)
@constraint(model, Oⁿˡ[3] -Pˢᴳₘₐₓ[3]*μᵐᵃˣˢᴳ³[T] +Pˢᴳₘᵢₙ[3]*μᵐⁱⁿˢᴳ³[T] +Kˢᵗ[3]*σˢᵗˢᴳ³[T] -Kˢʰ[3]*σˢʰˢᴳ³[T]+ ψᵐᵃˣˢᴳ³[T] >=0)
@constraint(model, Oⁿˡ[4] -Pˢᴳₘₐₓ[4]*μᵐᵃˣˢᴳ⁴[T] +Pˢᴳₘᵢₙ[4]*μᵐⁱⁿˢᴳ⁴[T] +Kˢᵗ[4]*σˢᵗˢᴳ⁴[T] -Kˢʰ[4]*σˢʰˢᴳ⁴[T]+ ψᵐᵃˣˢᴳ⁴[T] >=0)
for t in 1:T-1                                                                                                                                                       # dual constraints for UC, when t<=T-1
    @constraint(model, Oⁿˡ[1] -Pˢᴳₘₐₓ[1]*μᵐᵃˣˢᴳ¹[t]+ Pˢᴳₘᵢₙ[1]*μᵐⁱⁿˢᴳ¹[t]+ Kˢᵗ[1]*(σˢᵗˢᴳ¹[t]-σˢᵗˢᴳ¹[t+1])+ Kˢʰ[1]*(σˢʰˢᴳ¹[t+1]-σˢʰˢᴳ¹[t])+ ψᵐᵃˣˢᴳ¹[t] >=0)  
    @constraint(model, Oⁿˡ[2] -Pˢᴳₘₐₓ[2]*μᵐᵃˣˢᴳ²[t]+ Pˢᴳₘᵢₙ[2]*μᵐⁱⁿˢᴳ²[t]+ Kˢᵗ[2]*(σˢᵗˢᴳ²[t]-σˢᵗˢᴳ²[t+1])+ Kˢʰ[2]*(σˢʰˢᴳ²[t+1]-σˢʰˢᴳ²[t])+ ψᵐᵃˣˢᴳ²[t] >=0)
    @constraint(model, Oⁿˡ[3] -Pˢᴳₘₐₓ[3]*μᵐᵃˣˢᴳ³[t]+ Pˢᴳₘᵢₙ[3]*μᵐⁱⁿˢᴳ³[t]+ Kˢᵗ[3]*(σˢᵗˢᴳ³[t]-σˢᵗˢᴳ³[t+1])+ Kˢʰ[3]*(σˢʰˢᴳ³[t+1]-σˢʰˢᴳ³[t])+ ψᵐᵃˣˢᴳ³[t] >=0)
    @constraint(model, Oⁿˡ[4] -Pˢᴳₘₐₓ[4]*μᵐᵃˣˢᴳ⁴[t]+ Pˢᴳₘᵢₙ[4]*μᵐⁱⁿˢᴳ⁴[t]+ Kˢᵗ[4]*(σˢᵗˢᴳ⁴[t]-σˢᵗˢᴳ⁴[t+1])+ Kˢʰ[4]*(σˢʰˢᴳ⁴[t+1]-σˢʰˢᴳ⁴[t])+ ψᵐᵃˣˢᴳ⁴[t] >=0)
end

for t in 1:T                                                                                                  # dual constraints for generation, when t<=T-1, assume SGs in bus 2 are strategic
    @constraint(model, Oᵐ[1] -λᴱ[t] +μᵐᵃˣˢᴳ¹[t] -μᵐⁱⁿˢᴳ¹[t] >=0)                
    @constraint(model, Oᵐ[2]*b[t] -λᴱ[t] +μᵐᵃˣˢᴳ²[t] -μᵐⁱⁿˢᴳ²[t] >=0)                
    @constraint(model, Oᵐ[3] -λᴱ[t] +μᵐᵃˣˢᴳ³[t] -μᵐⁱⁿˢᴳ³[t] >=0)
    @constraint(model, Oᵐ[4] -λᴱ[t] +μᵐᵃˣˢᴳ⁴[t] -μᵐⁱⁿˢᴳ⁴[t] >=0)                
end
                                                                                          
@constraint(model,σˢᵗˢᴳ¹.<=1)                                 # dual constraints for on/off costs
@constraint(model,σˢᵗˢᴳ².<=1)
@constraint(model,σˢᵗˢᴳ³.<=1)
@constraint(model,σˢᵗˢᴳ⁴.<=1)



#-------Define Objective Functions 
@variable(model, Mc_1[1:T])   # auxiliary variable for revenue from market clearing
for t in 1:T
    @constraint(model, Mc_1[t]>= λᴱ[t]*Pˢᴳₘᵢₙ[2])
    @constraint(model, Mc_1[t]>= λᴱ[t]*Pˢᴳₘₐₓ[2] +bᵐᵃˣ*Oᵐ[2]*Pˢᴳ²[t] -bᵐᵃˣ*Oᵐ[2]*Pˢᴳₘₐₓ[2])
    @constraint(model, Mc_1[t]<= λᴱ[t]*Pˢᴳₘᵢₙ[2] +bᵐᵃˣ*Oᵐ[2]*Pˢᴳ²[t] -bᵐᵃˣ*Oᵐ[2]*Pˢᴳₘᵢₙ[2])
    @constraint(model, Mc_1[t]<= λᴱ[t]*Pˢᴳₘₐₓ[2] )
end
revenue_energy_UL=sum( Mc_1 )                      # revenue from market clearing
cost_nl_UL=sum(Oⁿˡ[2].*yˢᴳ² )                  # no-load cost of strategic SGs
cost_gene_UL=sum(Oᵐ[2].*Pˢᴳ²)                   # generation cost of strategic SGs 
cost_onoff_UL=sum(Cᵁ²)+sum(Cᴰ²)               # on/off cost of strategic SGs
obj_UL=revenue_energy_UL  -cost_nl_UL -cost_gene_UL -cost_onoff_UL                                            # objective function of UL

@variable(model, Mc_2[1:T])   # auxiliary variable 
for t in 1:T
    @constraint(model, Mc_2[t]>= b[t]*Pˢᴳₘᵢₙ[2] +Pˢᴳ²[t] -Pˢᴳₘᵢₙ[2])
    @constraint(model, Mc_2[t]>= b[t]*Pˢᴳₘₐₓ[2] +bᵐᵃˣ*Pˢᴳ²[t] -bᵐᵃˣ*Pˢᴳₘₐₓ[2])
    @constraint(model, Mc_2[t]<= b[t]*Pˢᴳₘᵢₙ[2] +bᵐᵃˣ*Pˢᴳ²[t] -bᵐᵃˣ*Pˢᴳₘᵢₙ[2])
    @constraint(model, Mc_2[t]<= b[t]*Pˢᴳₘₐₓ[2] +Pˢᴳ²[t] -Pˢᴳₘₐₓ[2])
end

cost_onoff_LL=sum(Cᵁ¹)+sum(Cᴰ¹)+sum(Cᵁ²)+sum(Cᴰ²)+sum(Cᵁ³)+sum(Cᴰ³)+sum(Cᵁ⁴)+sum(Cᴰ⁴)         
cost_nl_LL=sum(Oⁿˡ[1].*(yˢᴳ¹))+sum(Oⁿˡ[2].*(yˢᴳ²))+sum(Oⁿˡ[3].*(yˢᴳ³))+sum(Oⁿˡ[4].*(yˢᴳ⁴))   
cost_gene_LL=sum(Oᵐ[1]*Pˢᴳ¹) +sum(Oᵐ[2].*Mc_2) +sum(Oᵐ[3].*Pˢᴳ³) +sum(Oᵐ[4].*Pˢᴳ⁴)  

obj_LL=cost_onoff_LL+cost_nl_LL+cost_gene_LL 

@variable(model, obj_DLL_1[1:T])
for t in 1:T
    @constraint(model, obj_DLL_1[t]== Load_total[t]*λᴱ[t] -ψᵐᵃˣˢᴳ¹[t] -ψᵐᵃˣˢᴳ²[t] -ψᵐᵃˣˢᴳ³[t] -ψᵐᵃˣˢᴳ⁴[t] )
end

obj_DLL_2=  (yˢᴳ₀[1]*Kˢʰ[1]*σˢʰˢᴳ¹[1] +yˢᴳ₀[2]*Kˢʰ[2]*σˢʰˢᴳ²[1] +yˢᴳ₀[3]*Kˢʰ[3]*σˢʰˢᴳ³[1] +yˢᴳ₀[4]*Kˢʰ[4]*σˢʰˢᴳ⁴[1] ) 
            -(yˢᴳ₀[1]*Kˢᵗ[1]*σˢᵗˢᴳ¹[1] +yˢᴳ₀[2]*Kˢᵗ[2]*σˢᵗˢᴳ²[1] +yˢᴳ₀[3]*Kˢᵗ[3]*σˢᵗˢᴳ³[1] +yˢᴳ₀[4]*Kˢᵗ[4]*σˢᵗˢᴳ⁴[1] )  

obj_DLL=sum( obj_DLL_1 ) +sum(obj_DLL_2) 

W=2
  
@objective(model, Max, obj_UL-W*(obj_LL-obj_DLL))  # single-level objective function
#-------Solve and Output Results
#set_optimizer_attribute(model, "MIPGap", 0.01)
optimize!(model)

obj_UL=JuMP.value(obj_UL)
obj_LL=JuMP.value(obj_LL)
obj_DLL=JuMP.value(obj_DLL)
DG=obj_LL-obj_DLL

r_DG=DG/obj_LL*100

#-----------economic metrics for strategic one
λᴱ=JuMP.value.(λᴱ)
sum(λᴱ)/T

b=JuMP.value.(b)
sum(b)/T


yˢᴳ¹=JuMP.value.(yˢᴳ¹)
yˢᴳ²=JuMP.value.(yˢᴳ²)
yˢᴳ³=JuMP.value.(yˢᴳ³)
yˢᴳ⁴=JuMP.value.(yˢᴳ⁴)

Pˢᴳ¹=JuMP.value.(Pˢᴳ¹)
Pˢᴳ²=JuMP.value.(Pˢᴳ²)
Pˢᴳ³=JuMP.value.(Pˢᴳ³)
Pˢᴳ⁴=JuMP.value.(Pˢᴳ⁴)


non_RC_1=zeros(T)
non_RC_2=zeros(T)
non_RC_3=zeros(T)
non_RC_4=zeros(T)

RC_1=zeros(T)
RC_2=zeros(T)
RC_3=zeros(T)
RC_4=zeros(T)

for t in 1:T
    non_RC_1[t]=(Pˢᴳₘₐₓ[1] -Pˢᴳ¹[t])*(1 -yˢᴳ¹[t])
    non_RC_2[t]=(Pˢᴳₘₐₓ[2] -Pˢᴳ²[t])*(1 -yˢᴳ²[t])
    non_RC_3[t]=(Pˢᴳₘₐₓ[3] -Pˢᴳ³[t])*(1 -yˢᴳ³[t])
    non_RC_4[t]=(Pˢᴳₘₐₓ[4] -Pˢᴳ⁴[t])*(1 -yˢᴳ⁴[t])
    RC_1[t]=(Pˢᴳₘₐₓ[1] -Pˢᴳ¹[t])*yˢᴳ¹[t]
    RC_2[t]=(Pˢᴳₘₐₓ[2] -Pˢᴳ²[t])*yˢᴳ²[t]
    RC_3[t]=(Pˢᴳₘₐₓ[3] -Pˢᴳ³[t])*yˢᴳ³[t]
    RC_4[t]=(Pˢᴳₘₐₓ[4] -Pˢᴳ⁴[t])*yˢᴳ⁴[t]
end

plot(non_RC_1)
plot!(non_RC_2)
plot!(non_RC_3)
plot!(non_RC_4)

plot(RC_1)
plot!(RC_2)
plot!(RC_3)
plot!(RC_4)


