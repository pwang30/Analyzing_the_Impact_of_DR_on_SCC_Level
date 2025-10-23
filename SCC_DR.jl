# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# 19.Oct.2025

import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT
include("dataset_gene.jl")
include("offline_trainning.jl")
include("admittance_matrix_calculation.jl") 
# SGs, buses:2,3,4,5,27,30    IBRs, buses:1,23,26



#-----------------------------------Define Parameters for Calculating SCL-----------------------------------
I_IBG=1      # pre-defined SCL contribution from IBG
Iₗᵢₘ= 5       # SCL limit
β=0.95       # percentage of nominal voltage, range from 0.95-1.1
v_n=1        # nominal voltage
v=0.1        # gap for classification 

I_SCC_all_buses_scenarios, matrix_ω =dataset_gene(I_IBG, β,v_n)                                                            # data set generation                      
K_g, K_c, K_m, N_type_1, N_type_2, err_type_1, err_type_2= offline_trainning(I_SCC_all_buses_scenarios, matrix_ω, Iₗᵢₘ, v)  # offline_trainning

plot(I_SCC_appor[500:1500])
plot!(I_SCC_all_buses_scenarios[500:1500,1], seriestype=:scatter)

matwrite("IBG_bus1_real.mat", Dict("IBG_bus1_real" => I_SCC_all_buses_scenarios[500:1500,1]))

energy_price_restricted=[8.72 8.72 8.72 8.72 8.72 7.1 8.72 10.47 10.47 11.49 11.49 11.49 11.49 11.49 11.49 11.49 15.02 15.02 11.49 11.49 10.47 8.72 8.72 8.72 ]


#-----------------------------------Define Parameters for Optimization-----------------------------------  
      
Load_total=[18.42,17.95,18.29,18.51,18.13,17.88,19.46,21.97,23.17,23.87,
23.91,23.77,23.80,23.82,24.23,23.79,26.01,26.91,25.26,23.69,22.12,20.04,18.17,18.01]*10^3/3.5   # (MW)
T=length(Load_total)

IBG₁=250*ones(T,1)   
IBG₂₃=250*ones(T,1)
IBG₂₆=250*ones(T,1)


Pˢᴳₘₐₓ=[6.584, 5.760, 3.781, 3.335, 3.252, 2.880]*10^3/5            # Max generation of SGs                    SGs, buses:2,3,4,5,27,30
Pˢᴳₘᵢₙ=[3.292, 2.880, 1.512, 0.667, 0.650, 0.288]*10^3/5            #  Min generation of SGs                   SGs, buses:2,3,4,5,27,30
#Rₘₐₓ=[3.3, 2.9, 1.6, 1.334, 1.951, 1.728]*10^3/5                    #  Ramp limits of SGs                      SGs, buses:2,3,4,5,27,30
Kˢᵗ=[200, 125, 92.5, 72, 55, 31]*10^3 /10                           #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
Kˢʰ=[50, 28.5, 18.5, 14.4, 12, 10]*10^3  /10                        #  Shutdown cost of SGs                    SGs, buses:2,3,4,5,27,30
Oᵐ₁=[6.20, 7.10, 10.47, 12.28, 13.53, 15.36]                        #  Marginal generation cost of SG 1  in    SGs, buses:2,3,4,5,27,30
Oᵐ₂=[7.07, 8.72, 11.49, 12.84, 14.60, 15.02]                        #  Marginal generation cost of SG 2  in    SGs, buses:2,3,4,5,27,30
Oⁿˡ=[17.431, 15.005, 13.755, 10.930, 9.900, 8.570]*10^3 /10         #  No-load cost of SGs                     SGs, buses:2,3,4,5,27,30
P_g₀=[5.268 4.608 3.025 2.668 2.602 0]*10^3/5                       #  Initial generation (t=0) of SGs         SGs, buses:2,3,4,5,27,30
yˢᴳ₀=[1 1 1 1 1 0] 



#-----------------------------------Define Model-----------------------------------
model= Model()

@variable(model, Pˢᴳ²_1[1:T])        
@variable(model, Pˢᴳ²_2[1:T])  
@variable(model, Pˢᴳ³_1[1:T])     
@variable(model, Pˢᴳ³_2[1:T])                
@variable(model, Pˢᴳ⁴_1[1:T])   
@variable(model, Pˢᴳ⁴_2[1:T])             
@variable(model, Pˢᴳ⁵_1[1:T])    
@variable(model, Pˢᴳ⁵_2[1:T])            
@variable(model, Pˢᴳ²⁷_1[1:T]) 
@variable(model, Pˢᴳ²⁷_2[1:T])                
@variable(model, Pˢᴳ³⁰_1[1:T])                
@variable(model, Pˢᴳ³⁰_2[1:T])  

@variable(model, Pᴵᴮᴳ¹[1:T]>=0)         
@variable(model, Pᴵᴮᴳ²³[1:T]>=0)              
@variable(model, Pᴵᴮᴳ²⁶[1:T]>=0)  

@variable(model, yˢᴳ²_1[1:T], Bin)     
@variable(model, yˢᴳ²_2[1:T], Bin)          
@variable(model, yˢᴳ³_1[1:T], Bin)  
@variable(model, yˢᴳ³_2[1:T], Bin)          
@variable(model, yˢᴳ⁴_1[1:T], Bin)      
@variable(model, yˢᴳ⁴_2[1:T], Bin)      
@variable(model, yˢᴳ⁵_1[1:T], Bin)      
@variable(model, yˢᴳ⁵_2[1:T], Bin)      
@variable(model, yˢᴳ²⁷_1[1:T], Bin)  
@variable(model, yˢᴳ²⁷_2[1:T], Bin)          
@variable(model, yˢᴳ³⁰_1[1:T], Bin)   
@variable(model, yˢᴳ³⁰_2[1:T], Bin)                       


@variable(model, Cᵁ²_1[1:T]>=0)               
@variable(model, Cᵁ²_2[1:T]>=0)                 
@variable(model, Cᴰ²_1[1:T]>=0)                 
@variable(model, Cᴰ²_2[1:T]>=0)                 
@variable(model, Cᵁ³_1[1:T]>=0)    
@variable(model, Cᵁ³_2[1:T]>=0)             
@variable(model, Cᴰ³_1[1:T]>=0)  
@variable(model, Cᴰ³_2[1:T]>=0)               
@variable(model, Cᵁ⁴_1[1:T]>=0)     
@variable(model, Cᵁ⁴_2[1:T]>=0)            
@variable(model, Cᴰ⁴_1[1:T]>=0)   
@variable(model, Cᴰ⁴_2[1:T]>=0)              
@variable(model, Cᵁ⁵_1[1:T]>=0)  
@variable(model, Cᵁ⁵_2[1:T]>=0)               
@variable(model, Cᴰ⁵_1[1:T]>=0) 
@variable(model, Cᴰ⁵_2[1:T]>=0)                
@variable(model, Cᵁ²⁷_1[1:T]>=0)    
@variable(model, Cᵁ²⁷_2[1:T]>=0)             
@variable(model, Cᴰ²⁷_1[1:T]>=0)    
@variable(model, Cᴰ²⁷_2[1:T]>=0)             
@variable(model, Cᵁ³⁰_1[1:T]>=0)      
@variable(model, Cᵁ³⁰_2[1:T]>=0)            
@variable(model, Cᴰ³⁰_1[1:T]>=0)    
@variable(model, Cᴰ³⁰_2[1:T]>=0)            

@constraint(model, Pˢᴳ²_1.<=yˢᴳ²_1*Pˢᴳₘₐₓ[1])           # bounds for the output of SGs with UC , dual variables: μᵐⁱⁿₜ , μᵐᵃˣₜ
@constraint(model, yˢᴳ²_1*Pˢᴳₘᵢₙ[1].<=Pˢᴳ²_1)
@constraint(model, Pˢᴳ²_2.<=yˢᴳ²_2*Pˢᴳₘₐₓ[1])       
@constraint(model, yˢᴳ²_2*Pˢᴳₘᵢₙ[1].<=Pˢᴳ²_2)             
@constraint(model, Pˢᴳ³_1.<=yˢᴳ³_1*Pˢᴳₘₐₓ[2])       
@constraint(model, yˢᴳ³_1*Pˢᴳₘᵢₙ[2].<=Pˢᴳ³_1)       
@constraint(model, Pˢᴳ³_2.<=yˢᴳ³_2*Pˢᴳₘₐₓ[2])       
@constraint(model, yˢᴳ³_2*Pˢᴳₘᵢₙ[2].<=Pˢᴳ³_2)     
@constraint(model, Pˢᴳ⁴_1.<=yˢᴳ⁴_1*Pˢᴳₘₐₓ[3])       
@constraint(model, yˢᴳ⁴_1*Pˢᴳₘᵢₙ[3].<=Pˢᴳ⁴_1)
@constraint(model, Pˢᴳ⁴_2.<=yˢᴳ⁴_2*Pˢᴳₘₐₓ[3])       
@constraint(model, yˢᴳ⁴_2*Pˢᴳₘᵢₙ[3].<=Pˢᴳ⁴_2)     
@constraint(model, Pˢᴳ⁵_1.<=yˢᴳ⁵_1*Pˢᴳₘₐₓ[4])       
@constraint(model, yˢᴳ⁵_1*Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁵_1)
@constraint(model, Pˢᴳ⁵_2.<=yˢᴳ⁵_2*Pˢᴳₘₐₓ[4])       
@constraint(model, yˢᴳ⁵_2*Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁵_2)
@constraint(model, Pˢᴳ²⁷_1.<=yˢᴳ²⁷_1*Pˢᴳₘₐₓ[5])       
@constraint(model, yˢᴳ²⁷_1*Pˢᴳₘᵢₙ[5].<=Pˢᴳ²⁷_1)
@constraint(model, Pˢᴳ²⁷_2.<=yˢᴳ²⁷_2*Pˢᴳₘₐₓ[5])       
@constraint(model, yˢᴳ²⁷_2*Pˢᴳₘᵢₙ[5].<=Pˢᴳ²⁷_2)
@constraint(model, Pˢᴳ³⁰_1.<=yˢᴳ³⁰_1*Pˢᴳₘₐₓ[6])       
@constraint(model, yˢᴳ³⁰_1*Pˢᴳₘᵢₙ[6].<=Pˢᴳ³⁰_1)
@constraint(model, Pˢᴳ³⁰_2.<=yˢᴳ³⁰_2*Pˢᴳₘₐₓ[6])       
@constraint(model, yˢᴳ³⁰_2*Pˢᴳₘᵢₙ[6].<=Pˢᴳ³⁰_2)

@constraint(model, Cᵁ²_1[1]>=(yˢᴳ²_1[1]-yˢᴳ₀[1])*Kˢᵗ[1])        
@constraint(model, Cᴰ²_1[1]>=(yˢᴳ₀[1]-yˢᴳ²_1[1])*Kˢʰ[1])  
@constraint(model, Cᵁ²_2[1]>=(yˢᴳ²_2[1]-yˢᴳ₀[1])*Kˢᵗ[1])        
@constraint(model, Cᴰ²_2[1]>=(yˢᴳ₀[1]-yˢᴳ²_2[1])*Kˢʰ[1]) 
@constraint(model, Cᵁ³_1[1]>=(yˢᴳ³_1[1]-yˢᴳ₀[2])*Kˢᵗ[2])
@constraint(model, Cᴰ³_1[1]>=(yˢᴳ₀[2]-yˢᴳ³_1[1])*Kˢʰ[2])
@constraint(model, Cᵁ³_2[1]>=(yˢᴳ³_2[1]-yˢᴳ₀[2])*Kˢᵗ[2])
@constraint(model, Cᴰ³_2[1]>=(yˢᴳ₀[2]-yˢᴳ³_2[1])*Kˢʰ[2])
@constraint(model, Cᵁ⁴_1[1]>=(yˢᴳ⁴_1[1]-yˢᴳ₀[3])*Kˢᵗ[3])
@constraint(model, Cᴰ⁴_1[1]>=(yˢᴳ₀[3]-yˢᴳ⁴_1[1])*Kˢʰ[3])
@constraint(model, Cᵁ⁴_2[1]>=(yˢᴳ⁴_2[1]-yˢᴳ₀[3])*Kˢᵗ[3])
@constraint(model, Cᴰ⁴_2[1]>=(yˢᴳ₀[3]-yˢᴳ⁴_2[1])*Kˢʰ[3])
@constraint(model, Cᵁ⁵_1[1]>=(yˢᴳ⁵_1[1]-yˢᴳ₀[4])*Kˢᵗ[4])
@constraint(model, Cᴰ⁵_1[1]>=(yˢᴳ₀[4]-yˢᴳ⁵_1[1])*Kˢʰ[4])
@constraint(model, Cᵁ⁵_2[1]>=(yˢᴳ⁵_2[1]-yˢᴳ₀[4])*Kˢᵗ[4])
@constraint(model, Cᴰ⁵_2[1]>=(yˢᴳ₀[4]-yˢᴳ⁵_2[1])*Kˢʰ[4])
@constraint(model, Cᵁ²⁷_1[1]>=(yˢᴳ²⁷_1[1]-yˢᴳ₀[5])*Kˢᵗ[5])
@constraint(model, Cᴰ²⁷_1[1]>=(yˢᴳ₀[5]-yˢᴳ²⁷_1[1])*Kˢʰ[5])
@constraint(model, Cᵁ²⁷_2[1]>=(yˢᴳ²⁷_2[1]-yˢᴳ₀[5])*Kˢᵗ[5])
@constraint(model, Cᴰ²⁷_2[1]>=(yˢᴳ₀[5]-yˢᴳ²⁷_2[1])*Kˢʰ[5])
@constraint(model, Cᵁ³⁰_1[1]>=(yˢᴳ³⁰_1[1]-yˢᴳ₀[6])*Kˢᵗ[6])
@constraint(model, Cᴰ³⁰_1[1]>=(yˢᴳ₀[6]-yˢᴳ³⁰_1[1])*Kˢʰ[6])
@constraint(model, Cᵁ³⁰_2[1]>=(yˢᴳ³⁰_2[1]-yˢᴳ₀[6])*Kˢᵗ[6])
@constraint(model, Cᴰ³⁰_2[1]>=(yˢᴳ₀[6]-yˢᴳ³⁰_2[1])*Kˢʰ[6]) 
for t in 2:T
    @constraint(model, Cᵁ²_1[t]>=(yˢᴳ²_1[t]-yˢᴳ²_1[t-1])*Kˢᵗ[1])        
    @constraint(model, Cᴰ²_1[t]>=(yˢᴳ²_1[t-1]-yˢᴳ²_1[t])*Kˢʰ[1])  
    @constraint(model, Cᵁ²_2[t]>=(yˢᴳ²_2[t]-yˢᴳ²_2[t-1])*Kˢᵗ[1])        
    @constraint(model, Cᴰ²_2[t]>=(yˢᴳ²_2[t-1]-yˢᴳ²_2[t])*Kˢʰ[1]) 

    @constraint(model, Cᵁ³_1[t]>=(yˢᴳ³_1[t]-yˢᴳ³_1[t-1])*Kˢᵗ[2]) 
    @constraint(model, Cᴰ³_1[t]>=(yˢᴳ³_1[t-1]-yˢᴳ³_1[t])*Kˢʰ[2])
    @constraint(model, Cᵁ³_2[t]>=(yˢᴳ³_2[t]-yˢᴳ³_2[t-1])*Kˢᵗ[2])
    @constraint(model, Cᴰ³_2[t]>=(yˢᴳ³_2[t-1]-yˢᴳ³_2[t])*Kˢʰ[2])

    @constraint(model, Cᵁ⁴_1[t]>=(yˢᴳ⁴_1[t]-yˢᴳ⁴_1[t-1])*Kˢᵗ[3])
    @constraint(model, Cᴰ⁴_1[t]>=(yˢᴳ⁴_1[t-1]-yˢᴳ⁴_1[t])*Kˢʰ[3])
    @constraint(model, Cᵁ⁴_2[t]>=(yˢᴳ⁴_2[t]-yˢᴳ⁴_2[t-1])*Kˢᵗ[3])
    @constraint(model, Cᴰ⁴_2[t]>=(yˢᴳ⁴_2[t-1]-yˢᴳ⁴_2[t])*Kˢʰ[3])

    @constraint(model, Cᵁ⁵_1[t]>=(yˢᴳ⁵_1[t]-yˢᴳ⁵_1[t-1])*Kˢᵗ[4])
    @constraint(model, Cᴰ⁵_1[t]>=(yˢᴳ⁵_1[t-1]-yˢᴳ⁵_1[t])*Kˢʰ[4])
    @constraint(model, Cᵁ⁵_2[t]>=(yˢᴳ⁵_2[t]-yˢᴳ⁵_2[t-1])*Kˢᵗ[4])
    @constraint(model, Cᴰ⁵_2[t]>=(yˢᴳ⁵_2[t-1]-yˢᴳ⁵_2[t])*Kˢʰ[4])

    @constraint(model, Cᵁ²⁷_1[t]>=(yˢᴳ²⁷_1[t]-yˢᴳ²⁷_1[t-1])*Kˢᵗ[5])
    @constraint(model, Cᴰ²⁷_1[t]>=(yˢᴳ²⁷_1[t-1]-yˢᴳ²⁷_1[t])*Kˢʰ[5])
    @constraint(model, Cᵁ²⁷_2[t]>=(yˢᴳ²⁷_2[t]-yˢᴳ²⁷_2[t-1])*Kˢᵗ[5])
    @constraint(model, Cᴰ²⁷_2[t]>=(yˢᴳ²⁷_2[t-1]-yˢᴳ²⁷_2[t])*Kˢʰ[5])

    @constraint(model, Cᵁ³⁰_1[t]>=(yˢᴳ³⁰_1[t]-yˢᴳ³⁰_1[t-1])*Kˢᵗ[6])
    @constraint(model, Cᴰ³⁰_1[t]>=(yˢᴳ³⁰_1[t-1]-yˢᴳ³⁰_1[t])*Kˢʰ[6])
    @constraint(model, Cᵁ³⁰_2[t]>=(yˢᴳ³⁰_2[t]-yˢᴳ³⁰_2[t-1])*Kˢᵗ[6])
    @constraint(model, Cᴰ³⁰_2[t]>=(yˢᴳ³⁰_2[t-1]-yˢᴳ³⁰_2[t])*Kˢʰ[6])
end
    
for t in 1:T
    @constraint(model, Pᴵᴮᴳ¹[t] <= IBG₁[t])        
    @constraint(model, Pᴵᴮᴳ²³[t]<= IBG₂₃[t])       
    @constraint(model, Pᴵᴮᴳ²⁶[t]<= IBG₂₆[t])                  
end


@variable(model, I_scc[1:30, 1:T])  
for k in 1:30
for t in 1:T  
    @constraint(model,                                          
   I_scc[k,t] == K_g[1,k]*yˢᴳ²_1[t]+ K_g[2,k]*yˢᴳ²_2[t]+ 
    K_g[3,k]*yˢᴳ³_1[t]+ K_g[4,k]*yˢᴳ³_2[t]+ 
    K_g[5,k]*yˢᴳ⁴_1[t]+ K_g[6,k]*yˢᴳ⁴_2[t]+
    K_g[7,k]*yˢᴳ⁵_1[t]+ K_g[8,k]*yˢᴳ⁵_2[t]+
    K_g[9,k]*yˢᴳ²⁷_1[t]+ K_g[10,k]*yˢᴳ²⁷_2[t]+ 
    K_g[11,k]*yˢᴳ³⁰_1[t]+ K_g[12,k]*yˢᴳ³⁰_2[t]+

    K_c[1,k]+ K_c[2,k]+ K_c[3,k]+
    
    K_m[1,k]*yˢᴳ²_1[t]*yˢᴳ²_2[t]+ K_m[2,k]*yˢᴳ²_1[t]*yˢᴳ³_1[t]+ K_m[3,k]*yˢᴳ²_1[t]*yˢᴳ³_2[t]+ K_m[4,k]*yˢᴳ²_1[t]*yˢᴳ⁴_1[t]+
    K_m[5,k]*yˢᴳ²_1[t]*yˢᴳ⁴_2[t]+ K_m[6,k]*yˢᴳ²_1[t]*yˢᴳ⁵_1[t]+ K_m[7,k]*yˢᴳ²_1[t]*yˢᴳ⁵_2[t]+ K_m[8,k]*yˢᴳ²_1[t]*yˢᴳ²⁷_1[t]+
    K_m[9,k]*yˢᴳ²_1[t]*yˢᴳ²⁷_2[t]+ K_m[10,k]*yˢᴳ²_1[t]*yˢᴳ³⁰_1[t]+ K_m[11,k]*yˢᴳ²_1[t]*yˢᴳ³⁰_2[t]+

    K_m[12,k]*yˢᴳ²_2[t]*yˢᴳ³_1[t]+ K_m[13,k]*yˢᴳ²_2[t]*yˢᴳ³_2[t]+ K_m[14,k]*yˢᴳ²_2[t]*yˢᴳ⁴_1[t]+
    K_m[15,k]*yˢᴳ²_2[t]*yˢᴳ⁴_2[t]+ K_m[16,k]*yˢᴳ²_2[t]*yˢᴳ⁵_1[t]+ K_m[17,k]*yˢᴳ²_2[t]*yˢᴳ⁵_2[t]+ K_m[18,k]*yˢᴳ²_2[t]*yˢᴳ²⁷_1[t]+
    K_m[19,k]*yˢᴳ²_2[t]*yˢᴳ²⁷_2[t]+ K_m[20,k]*yˢᴳ²_2[t]*yˢᴳ³⁰_1[t]+ K_m[21,k]*yˢᴳ²_2[t]*yˢᴳ³⁰_2[t]+

    K_m[22,k]*yˢᴳ³_1[t]*yˢᴳ³_2[t]+ K_m[23,k]*yˢᴳ³_1[t]*yˢᴳ⁴_1[t]+
    K_m[24,k]*yˢᴳ³_1[t]*yˢᴳ⁴_2[t]+ K_m[25,k]*yˢᴳ³_1[t]*yˢᴳ⁵_1[t]+ K_m[26,k]*yˢᴳ³_1[t]*yˢᴳ⁵_2[t]+ K_m[27,k]*yˢᴳ³_1[t]*yˢᴳ²⁷_1[t]+
    K_m[28,k]*yˢᴳ³_1[t]*yˢᴳ²⁷_2[t]+ K_m[29,k]*yˢᴳ³_1[t]*yˢᴳ³⁰_1[t]+ K_m[30,k]*yˢᴳ³_1[t]*yˢᴳ³⁰_2[t]+

    K_m[31,k]*yˢᴳ³_2[t]*yˢᴳ⁴_1[t]+
    K_m[32,k]*yˢᴳ³_2[t]*yˢᴳ⁴_2[t]+ K_m[33,k]*yˢᴳ³_2[t]*yˢᴳ⁵_1[t]+ K_m[34,k]*yˢᴳ³_2[t]*yˢᴳ⁵_2[t]+ K_m[35,k]*yˢᴳ³_2[t]*yˢᴳ²⁷_1[t]+
    K_m[36,k]*yˢᴳ³_2[t]*yˢᴳ²⁷_2[t]+ K_m[37,k]*yˢᴳ³_2[t]*yˢᴳ³⁰_1[t]+ K_m[38,k]*yˢᴳ³_2[t]*yˢᴳ³⁰_2[t]+

    K_m[39,k]*yˢᴳ⁴_1[t]*yˢᴳ⁴_2[t]+ K_m[40,k]*yˢᴳ⁴_1[t]*yˢᴳ⁵_1[t]+ K_m[41,k]*yˢᴳ⁴_1[t]*yˢᴳ⁵_2[t]+ K_m[42,k]*yˢᴳ⁴_1[t]*yˢᴳ²⁷_1[t]+
    K_m[43,k]*yˢᴳ⁴_1[t]*yˢᴳ²⁷_2[t]+ K_m[44,k]*yˢᴳ⁴_1[t]*yˢᴳ³⁰_1[t]+ K_m[45,k]*yˢᴳ⁴_1[t]*yˢᴳ³⁰_2[t]+

    K_m[46,k]*yˢᴳ⁴_2[t]*yˢᴳ⁵_1[t]+ K_m[47,k]*yˢᴳ⁴_2[t]*yˢᴳ⁵_2[t]+ K_m[48,k]*yˢᴳ⁴_2[t]*yˢᴳ²⁷_1[t]+
    K_m[49,k]*yˢᴳ⁴_2[t]*yˢᴳ²⁷_2[t]+ K_m[50,k]*yˢᴳ⁴_2[t]*yˢᴳ³⁰_1[t]+ K_m[51,k]*yˢᴳ⁴_2[t]*yˢᴳ³⁰_2[t]+

    K_m[52,k]*yˢᴳ⁵_1[t]*yˢᴳ⁵_2[t]+ K_m[53,k]*yˢᴳ⁵_1[t]*yˢᴳ²⁷_1[t]+
    K_m[54,k]*yˢᴳ⁵_1[t]*yˢᴳ²⁷_2[t]+ K_m[55,k]*yˢᴳ⁵_1[t]*yˢᴳ³⁰_1[t]+ K_m[56,k]*yˢᴳ⁵_1[t]*yˢᴳ³⁰_2[t]+

    K_m[57,k]*yˢᴳ⁵_2[t]*yˢᴳ²⁷_1[t]+
    K_m[58,k]*yˢᴳ⁵_2[t]*yˢᴳ²⁷_2[t]+ K_m[59,k]*yˢᴳ⁵_2[t]*yˢᴳ³⁰_1[t]+ K_m[60,k]*yˢᴳ⁵_2[t]*yˢᴳ³⁰_2[t]+

    K_m[61,k]*yˢᴳ²⁷_1[t]*yˢᴳ²⁷_2[t]+ K_m[62,k]*yˢᴳ²⁷_1[t]*yˢᴳ³⁰_1[t]+ K_m[63,k]*yˢᴳ²⁷_1[t]*yˢᴳ³⁰_2[t]+

    K_m[64,k]*yˢᴳ²⁷_2[t]*yˢᴳ³⁰_1[t]+ K_m[65,k]*yˢᴳ²⁷_2[t]*yˢᴳ³⁰_2[t]+

    K_m[66,k]*yˢᴳ³⁰_1[t]*yˢᴳ³⁰_2[t] )

    @constraint(model,  I_scc[k,t] >= Iₗᵢₘ )

end
end


#--------DR modelling

xs = energy_price_restricted

kil = [50, 70, 100]         # Interrupted Load Compensation Fee
cil = [0.1, 0.08, 0.05]      # Maximum Interrupted Load Ratio
sil=[30, 20]  

@variable(model, 0 <= pil[1:3, 1:T])  # Interrupted Load
@variable(model, 0 <= shiftl[1:T])    # Transferred-out Load
@variable(model, 0 <= shiftq[1:T])    # Transferred-in Load
@variable(model, ushiftl[1:T], Bin)   
@variable(model, ushiftq[1:T], Bin)
@variable(model, pmgb[1:T])           # The User's Final Electricity Consumption
@variable(model, cshift[1:T])         # Transferred Load Compensation Fee

for t in 1:T
    @constraint(model, Pˢᴳ²_1[t]+Pˢᴳ²_2[t]+Pˢᴳ³_1[t]+Pˢᴳ³_2[t]+Pˢᴳ⁴_1[t]+Pˢᴳ⁴_2[t]+Pˢᴳ⁵_1[t]+Pˢᴳ⁵_2[t]+Pˢᴳ²⁷_1[t]+Pˢᴳ²⁷_2[t]+Pˢᴳ³⁰_1[t]+Pˢᴳ³⁰_2[t]+
                   Pᴵᴮᴳ¹[t]+Pᴵᴮᴳ²³[t]+Pᴵᴮᴳ²⁶[t] == pmgb[t])     # power balance 
end

pload=Load_total   # (MW)
for m in 1:3, t in 1:T   # Hierarchical Interrupted Load Constraints
    @constraint(model, pil[m,t] <= cil[m]*pload[t])
end

for m in 1:3, t in 2:T   # Continuous Interrupted Load Constraints
    @constraint(model, pil[m,t] + pil[m,t-1] <= 0.2*pload[t])
end

for t in 1:T               # Transferable Load Constraints
    @constraint(model, ushiftl[t] + ushiftq[t] <= 1)
    @constraint(model, shiftl[t] <= ushiftl[t]*pload[t]*0.12)
    @constraint(model, shiftq[t] <= ushiftq[t]*pload[t]*0.12)
end

@constraint(model, sum(shiftl) == sum(shiftq))  # Transfer Load Conservation

for t in 1:T               # Power Balance Constraints
    @constraint(model, pload[t] - shiftl[t] + shiftq[t] - sum(pil[:,t]) == pmgb[t])
end

for t in 1:T           # Calculate the Compensation Fee for Transferred Load
    @constraint(model, cshift[t] == sil[1]*shiftl[t] + sil[2]*shiftq[t])
end

@variable(model, Cost_of_users)   # Cost of electricity users
@constraint(model,  Cost_of_users == sum(xs.*pmgb ) - kil[1]*(sum(pil[1,:])) - kil[2]*(sum(pil[2,:])) - kil[3]*(sum(pil[3,:])) - sum(cshift))


cost_onoff_Primal=sum(Cᵁ²_1)+sum(Cᴰ²_1)+sum(Cᵁ³_1)+sum(Cᴰ³_1)+sum(Cᵁ⁴_1)+sum(Cᴰ⁴_1)+sum(Cᵁ⁵_1)+sum(Cᴰ⁵_1)+sum(Cᵁ²⁷_1)+sum(Cᴰ²⁷_1)+sum(Cᵁ³⁰_1)+sum(Cᴰ³⁰_1)  +sum(Cᵁ²_2)+sum(Cᴰ²_2)+sum(Cᵁ³_2)+sum(Cᴰ³_2)+sum(Cᵁ⁴_2)+sum(Cᴰ⁴_2)+sum(Cᵁ⁵_2)+sum(Cᴰ⁵_2)+sum(Cᵁ²⁷_2)+sum(Cᴰ²⁷_2)+sum(Cᵁ³⁰_2)+sum(Cᴰ³⁰_2)       
cost_nl_Primal=sum(Oⁿˡ[1].*(yˢᴳ²_1+yˢᴳ²_2))+sum(Oⁿˡ[2].*(yˢᴳ³_1+yˢᴳ³_2))+sum(Oⁿˡ[3].*(yˢᴳ⁴_1+yˢᴳ⁴_2))+sum(Oⁿˡ[4].*(yˢᴳ⁵_1+yˢᴳ⁵_2))+sum(Oⁿˡ[5].*(yˢᴳ²⁷_1+yˢᴳ²⁷_2))+sum(Oⁿˡ[6].*(yˢᴳ³⁰_1+yˢᴳ³⁰_2))    
cost_gene_Primal=sum(Oᵐ₁[1].*Pˢᴳ²_1+Oᵐ₂[1].*Pˢᴳ²_2 )+ sum(Oᵐ₁[2].*Pˢᴳ³_1+Oᵐ₂[2].*Pˢᴳ³_2 )+sum(Oᵐ₁[3].*Pˢᴳ⁴_1+Oᵐ₂[3].*Pˢᴳ⁴_2)+sum(Oᵐ₁[4].*Pˢᴳ⁵_1+Oᵐ₂[4].*Pˢᴳ⁵_2)+sum(Oᵐ₁[5].*Pˢᴳ²⁷_1+Oᵐ₂[5].*Pˢᴳ²⁷_2)+sum(Oᵐ₁[6].*Pˢᴳ³⁰_1+Oᵐ₂[6].*Pˢᴳ³⁰_2)   

obj_Primal=cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal


@objective(model, Min, obj_Primal + Cost_of_users)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
optimize!(model)

pil= value.(pil)

obj_Primal_value= value.(obj_Primal)
Cost_of_users_value= value.(Cost_of_users)
obj_Primal_value+Cost_of_users_value

matwrite("energy_price.mat", Dict("energy_price" => energy_price_restricted))


Cost_of_users_value= value.(Cost_of_users)

yˢᴳ²_1= value.(yˢᴳ²_1)
yˢᴳ²_2= value.(yˢᴳ²_2)
yˢᴳ³_1= value.(yˢᴳ³_1)
yˢᴳ³_2= value.(yˢᴳ³_2)
yˢᴳ⁴_1= value.(yˢᴳ⁴_1)
yˢᴳ⁴_2= value.(yˢᴳ⁴_2)
yˢᴳ⁵_1= value.(yˢᴳ⁵_1)
yˢᴳ⁵_2= value.(yˢᴳ⁵_2)
yˢᴳ²⁷_1= value.(yˢᴳ²⁷_1)
yˢᴳ²⁷_2= value.(yˢᴳ²⁷_2)
yˢᴳ³⁰_1= value.(yˢᴳ³⁰_1)
yˢᴳ³⁰_2= value.(yˢᴳ³⁰_2)


I_min=zeros(1,30)
I_scc=zeros(30,T)

for k in 1:30
for t in 1:T                                             # bounds for the SCL of buses  I_₂₆   
   I_scc[k,t]=K_g[1,k]*yˢᴳ²_1[t]+ K_g[2,k]*yˢᴳ²_2[t]+ 
    K_g[3,k]*yˢᴳ³_1[t]+ K_g[4,k]*yˢᴳ³_2[t]+ 
    K_g[5,k]*yˢᴳ⁴_1[t]+ K_g[6,k]*yˢᴳ⁴_2[t]+
    K_g[7,k]*yˢᴳ⁵_1[t]+ K_g[8,k]*yˢᴳ⁵_2[t]+
    K_g[9,k]*yˢᴳ²⁷_1[t]+ K_g[10,k]*yˢᴳ²⁷_2[t]+ 
    K_g[11,k]*yˢᴳ³⁰_1[t]+ K_g[12,k]*yˢᴳ³⁰_2[t]+

    K_c[1,k]+ K_c[2,k]+ K_c[3,k]+
    
    K_m[1,k]*yˢᴳ²_1[t]*yˢᴳ²_2[t]+ K_m[2,k]*yˢᴳ²_1[t]*yˢᴳ³_1[t]+ K_m[3,k]*yˢᴳ²_1[t]*yˢᴳ³_2[t]+ K_m[4,k]*yˢᴳ²_1[t]*yˢᴳ⁴_1[t]+
    K_m[5,k]*yˢᴳ²_1[t]*yˢᴳ⁴_2[t]+ K_m[6,k]*yˢᴳ²_1[t]*yˢᴳ⁵_1[t]+ K_m[7,k]*yˢᴳ²_1[t]*yˢᴳ⁵_2[t]+ K_m[8,k]*yˢᴳ²_1[t]*yˢᴳ²⁷_1[t]+
    K_m[9,k]*yˢᴳ²_1[t]*yˢᴳ²⁷_2[t]+ K_m[10,k]*yˢᴳ²_1[t]*yˢᴳ³⁰_1[t]+ K_m[11,k]*yˢᴳ²_1[t]*yˢᴳ³⁰_2[t]+

    K_m[12,k]*yˢᴳ²_2[t]*yˢᴳ³_1[t]+ K_m[13,k]*yˢᴳ²_2[t]*yˢᴳ³_2[t]+ K_m[14,k]*yˢᴳ²_2[t]*yˢᴳ⁴_1[t]+
    K_m[15,k]*yˢᴳ²_2[t]*yˢᴳ⁴_2[t]+ K_m[16,k]*yˢᴳ²_2[t]*yˢᴳ⁵_1[t]+ K_m[17,k]*yˢᴳ²_2[t]*yˢᴳ⁵_2[t]+ K_m[18,k]*yˢᴳ²_2[t]*yˢᴳ²⁷_1[t]+
    K_m[19,k]*yˢᴳ²_2[t]*yˢᴳ²⁷_2[t]+ K_m[20,k]*yˢᴳ²_2[t]*yˢᴳ³⁰_1[t]+ K_m[21,k]*yˢᴳ²_2[t]*yˢᴳ³⁰_2[t]+

    K_m[22,k]*yˢᴳ³_1[t]*yˢᴳ³_2[t]+ K_m[23,k]*yˢᴳ³_1[t]*yˢᴳ⁴_1[t]+
    K_m[24,k]*yˢᴳ³_1[t]*yˢᴳ⁴_2[t]+ K_m[25,k]*yˢᴳ³_1[t]*yˢᴳ⁵_1[t]+ K_m[26,k]*yˢᴳ³_1[t]*yˢᴳ⁵_2[t]+ K_m[27,k]*yˢᴳ³_1[t]*yˢᴳ²⁷_1[t]+
    K_m[28,k]*yˢᴳ³_1[t]*yˢᴳ²⁷_2[t]+ K_m[29,k]*yˢᴳ³_1[t]*yˢᴳ³⁰_1[t]+ K_m[30,k]*yˢᴳ³_1[t]*yˢᴳ³⁰_2[t]+

    K_m[31,k]*yˢᴳ³_2[t]*yˢᴳ⁴_1[t]+
    K_m[32,k]*yˢᴳ³_2[t]*yˢᴳ⁴_2[t]+ K_m[33,k]*yˢᴳ³_2[t]*yˢᴳ⁵_1[t]+ K_m[34,k]*yˢᴳ³_2[t]*yˢᴳ⁵_2[t]+ K_m[35,k]*yˢᴳ³_2[t]*yˢᴳ²⁷_1[t]+
    K_m[36,k]*yˢᴳ³_2[t]*yˢᴳ²⁷_2[t]+ K_m[37,k]*yˢᴳ³_2[t]*yˢᴳ³⁰_1[t]+ K_m[38,k]*yˢᴳ³_2[t]*yˢᴳ³⁰_2[t]+

    K_m[39,k]*yˢᴳ⁴_1[t]*yˢᴳ⁴_2[t]+ K_m[40,k]*yˢᴳ⁴_1[t]*yˢᴳ⁵_1[t]+ K_m[41,k]*yˢᴳ⁴_1[t]*yˢᴳ⁵_2[t]+ K_m[42,k]*yˢᴳ⁴_1[t]*yˢᴳ²⁷_1[t]+
    K_m[43,k]*yˢᴳ⁴_1[t]*yˢᴳ²⁷_2[t]+ K_m[44,k]*yˢᴳ⁴_1[t]*yˢᴳ³⁰_1[t]+ K_m[45,k]*yˢᴳ⁴_1[t]*yˢᴳ³⁰_2[t]+

    K_m[46,k]*yˢᴳ⁴_2[t]*yˢᴳ⁵_1[t]+ K_m[47,k]*yˢᴳ⁴_2[t]*yˢᴳ⁵_2[t]+ K_m[48,k]*yˢᴳ⁴_2[t]*yˢᴳ²⁷_1[t]+
    K_m[49,k]*yˢᴳ⁴_2[t]*yˢᴳ²⁷_2[t]+ K_m[50,k]*yˢᴳ⁴_2[t]*yˢᴳ³⁰_1[t]+ K_m[51,k]*yˢᴳ⁴_2[t]*yˢᴳ³⁰_2[t]+

    K_m[52,k]*yˢᴳ⁵_1[t]*yˢᴳ⁵_2[t]+ K_m[53,k]*yˢᴳ⁵_1[t]*yˢᴳ²⁷_1[t]+
    K_m[54,k]*yˢᴳ⁵_1[t]*yˢᴳ²⁷_2[t]+ K_m[55,k]*yˢᴳ⁵_1[t]*yˢᴳ³⁰_1[t]+ K_m[56,k]*yˢᴳ⁵_1[t]*yˢᴳ³⁰_2[t]+

    K_m[57,k]*yˢᴳ⁵_2[t]*yˢᴳ²⁷_1[t]+
    K_m[58,k]*yˢᴳ⁵_2[t]*yˢᴳ²⁷_2[t]+ K_m[59,k]*yˢᴳ⁵_2[t]*yˢᴳ³⁰_1[t]+ K_m[60,k]*yˢᴳ⁵_2[t]*yˢᴳ³⁰_2[t]+

    K_m[61,k]*yˢᴳ²⁷_1[t]*yˢᴳ²⁷_2[t]+ K_m[62,k]*yˢᴳ²⁷_1[t]*yˢᴳ³⁰_1[t]+ K_m[63,k]*yˢᴳ²⁷_1[t]*yˢᴳ³⁰_2[t]+

    K_m[64,k]*yˢᴳ²⁷_2[t]*yˢᴳ³⁰_1[t]+ K_m[65,k]*yˢᴳ²⁷_2[t]*yˢᴳ³⁰_2[t]+

    K_m[66,k]*yˢᴳ³⁰_1[t]*yˢᴳ³⁰_2[t]

end
I_min[k]=minimum(I_scc[k,:])
end


matwrite("I_min_no_DR.mat", Dict("I_min_no_DR" => I_min))





