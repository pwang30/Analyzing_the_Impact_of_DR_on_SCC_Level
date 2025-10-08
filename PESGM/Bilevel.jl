# PESGM 2026, Canada 
# 19.Jul.2025 

import Pkg
#Xusing JuMP,Gurobi, Bilevel, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT, CPLEX, Ipopt
using JuMP, BilevelJuMP,Plots,MAT, Ipopt


#----------------------------------Parameters----------------------------------#
Load_orig=[18.42,17.95,18.29,18.51,18.13,17.88,19.46,21.97,23.17,23.87,
23.91,23.77,23.80,23.82,24.23,23.79,26.01,26.91,25.26,23.69,22.12,20.04,18.17,18.01]*10   # (MW)
T=length(Load_orig)

Pˢᴳₘₐₓ=[6.584, 5.760, 3.781, 3.335, 3.252, 2.880]*15            
Pˢᴳₘᵢₙ=[3.292, 2.880, 1.512, 0.667, 0.650, 0.288]*15            
Rₘₐₓ=[1.317, 1.152, 1.512, 1.334, 1.951, 1.728]*15              
Oᵐ=[6.20, 7.10, 10.47, 12.28, 13.53, 15.36]                        
P_g₀=[5.268 4.608 3.025 2.668 2.602 0]*15                      

k_min=1
k_max=2

penalty=30
comp1=100
comp2=30


#-----------------------------------Define Bilevel Model-----------------------------------
model= BilevelModel()
 
@variable(Lower(model),  Pˢᴳ²[1:T])          
@variable(Lower(model),  Pˢᴳ³[1:T])     
@variable(Lower(model),  Pˢᴳ⁴[1:T])              
@variable(Lower(model),  Pˢᴳ⁵[1:T])            
@variable(Lower(model),  Pˢᴳ²⁷[1:T])               
@variable(Lower(model),  Pˢᴳ³⁰[1:T])                
                  

@constraint(Lower(model),  Pˢᴳ².<=Pˢᴳₘₐₓ[1])       
@constraint(Lower(model),  Pˢᴳₘᵢₙ[1].<=Pˢᴳ²)             
@constraint(Lower(model),  Pˢᴳ³.<=Pˢᴳₘₐₓ[2])       
@constraint(Lower(model),  Pˢᴳₘᵢₙ[2].<=Pˢᴳ³)           
@constraint(Lower(model),  Pˢᴳ⁴.<=Pˢᴳₘₐₓ[3])       
@constraint(Lower(model),  Pˢᴳₘᵢₙ[3].<=Pˢᴳ⁴)    
@constraint(Lower(model),  Pˢᴳ⁵.<=Pˢᴳₘₐₓ[4])       
@constraint(Lower(model),  Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁵)
@constraint(Lower(model),  Pˢᴳ²⁷.<=Pˢᴳₘₐₓ[5])       
@constraint(Lower(model),  Pˢᴳₘᵢₙ[5].<=Pˢᴳ²⁷)
@constraint(Lower(model),  Pˢᴳ³⁰.<=Pˢᴳₘₐₓ[6])       
@constraint(Lower(model),  Pˢᴳₘᵢₙ[6].<=Pˢᴳ³⁰)

@constraint(Lower(model),  Pˢᴳ²[1]-P_g₀[1]<=Rₘₐₓ[1])       
@constraint(Lower(model),  -Rₘₐₓ[1]<=Pˢᴳ²[1]-P_g₀[1])  
@constraint(Lower(model),  Pˢᴳ³[1]-P_g₀[2]<=Rₘₐₓ[2])        
@constraint(Lower(model),  -Rₘₐₓ[2]<=Pˢᴳ³[1]-P_g₀[2]) 
@constraint(Lower(model),  Pˢᴳ⁴[1]-P_g₀[3]<=Rₘₐₓ[3])
@constraint(Lower(model),  -Rₘₐₓ[3]<=Pˢᴳ⁴[1]-P_g₀[3])
@constraint(Lower(model),  Pˢᴳ⁵[1]-P_g₀[4]<=Rₘₐₓ[4])
@constraint(Lower(model),  -Rₘₐₓ[4]<=Pˢᴳ⁵[1]-P_g₀[4])
@constraint(Lower(model),  Pˢᴳ²⁷[1]-P_g₀[5]<=Rₘₐₓ[5])
@constraint(Lower(model),  -Rₘₐₓ[5]<=Pˢᴳ²⁷[1]-P_g₀[5])
@constraint(Lower(model),  Pˢᴳ³⁰[1]-P_g₀[6]<=Rₘₐₓ[6])
@constraint(Lower(model),  -Rₘₐₓ[6]<=Pˢᴳ³⁰[1]-P_g₀[6])

for t in 2:T                                                   
    @constraint(Lower(model),  Pˢᴳ²[t]-Pˢᴳ²[t-1]<=Rₘₐₓ[1])        
    @constraint(Lower(model),  -Rₘₐₓ[1]<=Pˢᴳ²[t]-Pˢᴳ²[t-1])  
    @constraint(Lower(model),  Pˢᴳ³[t]-Pˢᴳ³[t-1]<=Rₘₐₓ[2])
    @constraint(Lower(model),  -Rₘₐₓ[2]<=Pˢᴳ³[t]-Pˢᴳ³[t-1]) 
    @constraint(Lower(model),  Pˢᴳ⁴[t]-Pˢᴳ⁴[t-1]<=Rₘₐₓ[3])
    @constraint(Lower(model),  -Rₘₐₓ[3]<=Pˢᴳ⁴[t]-Pˢᴳ⁴[t-1])
    @constraint(Lower(model),  Pˢᴳ⁵[t]-Pˢᴳ⁵[t-1]<=Rₘₐₓ[4])
    @constraint(Lower(model),  -Rₘₐₓ[4]<=Pˢᴳ⁵[t]-Pˢᴳ⁵[t-1])
    @constraint(Lower(model),  Pˢᴳ²⁷[t]-Pˢᴳ²⁷[t-1]<=Rₘₐₓ[5])
    @constraint(Lower(model),  -Rₘₐₓ[5]<=Pˢᴳ²⁷[t]-Pˢᴳ²⁷[t-1])
    @constraint(Lower(model),  Pˢᴳ³⁰[t]-Pˢᴳ³⁰[t-1]<=Rₘₐₓ[6])
    @constraint(Lower(model),  -Rₘₐₓ[6]<=Pˢᴳ³⁰[t]-Pˢᴳ³⁰[t-1])
end


#---------------------Model of demand response---------------------

@variable(Lower(model),  P_D[1:T])
@constraint(Lower(model),  P_D >= 0)  
#@constraint(Lower(model),  P_D <= Load_orig) 
@constraint(Lower(model),  sum(P_D) >= 0.8*sum(Load_orig))

for t in 1:8
    @constraint(Lower(model),  P_D[t] >= Load_orig[t]+2 )
end
for t in 21:24
    @constraint(Lower(model),  P_D[t] >= Load_orig[t]+2 )
end
for t in 15:18
    @constraint(Lower(model),  P_D[t] <= Load_orig[t]-6 )
end

power_balance=Dict()
for t in 1:T
    power_balance[t]=@constraint(Lower(model), Pˢᴳ²[t]+Pˢᴳ³[t]+Pˢᴳ⁴[t]+Pˢᴳ⁵[t]+Pˢᴳ²⁷[t]+Pˢᴳ³⁰[t] == P_D[t])     
end

@variable(Upper(model), k_min <= k[1:T] <= k_max)
@variable(Upper(model), lambda_1<=15.36, DualOf(power_balance[1]))
@variable(Upper(model), lambda_2<=15.36, DualOf(power_balance[2]))
@variable(Upper(model), lambda_3<=15.36, DualOf(power_balance[3]))
@variable(Upper(model), lambda_4<=15.36, DualOf(power_balance[4]))
@variable(Upper(model), lambda_5<=15.36, DualOf(power_balance[5]))
@variable(Upper(model), lambda_6<=15.36, DualOf(power_balance[6]))
@variable(Upper(model), lambda_7<=15.36, DualOf(power_balance[7]))
@variable(Upper(model), lambda_8<=15.36, DualOf(power_balance[8]))
@variable(Upper(model), lambda_9<=15.36, DualOf(power_balance[9]))
@variable(Upper(model), lambda_10<=15.36, DualOf(power_balance[10]))
@variable(Upper(model), lambda_11<=15.36, DualOf(power_balance[11]))
@variable(Upper(model), lambda_12<=15.36, DualOf(power_balance[12]))
@variable(Upper(model), lambda_13<=15.36, DualOf(power_balance[13]))
@variable(Upper(model), lambda_14<=15.36, DualOf(power_balance[14]))
@variable(Upper(model), lambda_15<=15.36, DualOf(power_balance[15]))
@variable(Upper(model), lambda_16<=15.36, DualOf(power_balance[16]))
@variable(Upper(model), lambda_17<=15.36, DualOf(power_balance[17]))
@variable(Upper(model), lambda_18<=15.36, DualOf(power_balance[18]))
@variable(Upper(model), lambda_19<=15.36, DualOf(power_balance[19]))
@variable(Upper(model), lambda_20<=15.36, DualOf(power_balance[20]))
@variable(Upper(model), lambda_21<=15.36, DualOf(power_balance[21]))
@variable(Upper(model), lambda_22<=15.36, DualOf(power_balance[22]))
@variable(Upper(model), lambda_23<=15.36, DualOf(power_balance[23]))
@variable(Upper(model), lambda_24<=15.36, DualOf(power_balance[24]))

#-------Off-Peak  01:00-08:00 21:00-24:00.   Flat.  09:00-14:00    19:00-20:00.   Peak.  15:00-18:00




rev_comp_1=lambda_1*Pˢᴳ²[1] +lambda_2*Pˢᴳ²[2] +lambda_3*Pˢᴳ²[3] +lambda_4*Pˢᴳ²[4] +lambda_5*Pˢᴳ²[5] +lambda_6*Pˢᴳ²[6] +lambda_7*Pˢᴳ²[7] +lambda_8*Pˢᴳ²[8] +lambda_9*Pˢᴳ²[9] +lambda_10*Pˢᴳ²[10] + 
lambda_11*Pˢᴳ²[11] +lambda_12*Pˢᴳ²[12] +lambda_13*Pˢᴳ²[13] +lambda_14*Pˢᴳ²[14] +lambda_15*Pˢᴳ²[15] +lambda_16*Pˢᴳ²[16] +lambda_17*Pˢᴳ²[17] +lambda_18*Pˢᴳ²[18] +lambda_19*Pˢᴳ²[19] +lambda_20*Pˢᴳ²[20] +
lambda_21*Pˢᴳ²[21] +lambda_22*Pˢᴳ²[22] +lambda_23*Pˢᴳ²[23] +lambda_24*Pˢᴳ²[24]

@objective(Upper(model), Max, rev_comp_1 -(Oᵐ[1]*sum(Pˢᴳ²)) )


C_1=sum(k.*Pˢᴳ²*Oᵐ[1])
C_2=Oᵐ[2]*sum(Pˢᴳ³)
C_3=Oᵐ[3]*sum(Pˢᴳ⁴)
C_4=Oᵐ[4]*sum(Pˢᴳ⁵)
C_5=Oᵐ[5]*sum(Pˢᴳ²⁷)
C_6=Oᵐ[6]*sum(Pˢᴳ³⁰)

@objective(Lower(model), Min, C_1 + C_2 + C_3 + C_4 + C_5 + C_6  -(P_D[9]-Load_orig[9] +P_D[10]-Load_orig[10] +P_D[11]-Load_orig[11] 
                                +P_D[12] -Load_orig[12] +P_D[13] -Load_orig[13] +P_D[14] -Load_orig[14] +P_D[19] -Load_orig[19] 
                                +P_D[20] -Load_orig[20] )*penalty )

#@objective(Lower(model), Min, C_1 + C_2 + C_3 + C_4 + C_5 + C_6  -(P_D[1]-Load_orig[1] +P_D[2]-Load_orig[2] +P_D[3]-Load_orig[3] 
                                #+P_D[4] -Load_orig[4] +P_D[5] -Load_orig[5] +P_D[6] -Load_orig[6] +P_D[7] -Load_orig[7] 
                                #+P_D[8] -Load_orig[8] +P_D[21]-Load_orig[21] 
                                #+P_D[22]-Load_orig[22] +P_D[23]-Load_orig[23] +P_D[24]-Load_orig[24])*comp1 )



#---------------------------solve the bi-level model-------------------

BilevelJuMP.set_mode(model, BilevelJuMP.StrongDualityMode())
set_optimizer(model, Ipopt.Optimizer)
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


lambda_1=JuMP.value.(lambda_1)
lambda_2=JuMP.value.(lambda_2)
lambda_3=JuMP.value.(lambda_3)
lambda_4=JuMP.value.(lambda_4)
lambda_5=JuMP.value.(lambda_5)
lambda_6=JuMP.value.(lambda_6)
lambda_7=JuMP.value.(lambda_7)
lambda_8=JuMP.value.(lambda_8)
lambda_9=JuMP.value.(lambda_9)
lambda_10=JuMP.value.(lambda_10)
lambda_11=JuMP.value.(lambda_11)
lambda_12=JuMP.value.(lambda_12)
lambda_13=JuMP.value.(lambda_13)
lambda_14=JuMP.value.(lambda_14)
lambda_15=JuMP.value.(lambda_15)
lambda_16=JuMP.value.(lambda_16)
lambda_17=JuMP.value.(lambda_17)
lambda_18=JuMP.value.(lambda_18)
lambda_19=JuMP.value.(lambda_19)
lambda_20=JuMP.value.(lambda_20)
lambda_21=JuMP.value.(lambda_21)
lambda_22=JuMP.value.(lambda_22)
lambda_23=JuMP.value.(lambda_23)
lambda_24=JuMP.value.(lambda_24)
clearing_price=[lambda_1,lambda_2,lambda_3,lambda_4,lambda_5,lambda_6,lambda_7,lambda_8,lambda_9,lambda_10,
lambda_11,lambda_12,lambda_13,lambda_14,lambda_15,lambda_16,lambda_17,lambda_18,lambda_19,
lambda_20,lambda_21,lambda_22,lambda_23,lambda_24]

plot(k)
plot(clearing_price)

plot(Load_orig)
plot!(P_D)

matwrite("SCL_prices_bus30_disp_MAC.mat", Dict("SCL_prices_bus30_disp_MAC" => SCL_prices_bus30_MAC))