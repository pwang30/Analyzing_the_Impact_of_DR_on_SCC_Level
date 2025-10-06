# Analyzing the Role of DSO in Electricity Trading of VPPs via a Stackelberg Game Model
# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Note 1: This code is full, easy to understand and replicate. Every line has its corresponding notes for a better understanding.
# Note 2: For each model of VPP/Lower level, they have same codeing structure, the only difference is the parameters. 
#         So, once you know how to structur anyone of them, please dive into DSO/Upper level modeling.

# 19.Dec.2024

#----------------------Pkgs introduction-----------------------

import Pkg 
Pkg.add(BilevelJuMP)
Pkg.add(SCIP)
using JuMP,BilevelJuMP,SCIP,XLSX   

#----------------------define model----------------------------

model = BilevelModel(
   SCIP.Optimizer,
)

#----------------------relevant parameters input---------------

P_Wmax_1=[2,1.5,1.6,1.8,1.3,0.6,2.8,3.3,3.9,4,3.3,2.9,2.7,2,0.2,3.2,5.1,3.1,1.8,2,1.3,1,2,3.8]                    # Wind_VPP1
load_1=[2.2,1.8,3,6,5.8,5.2,5.6,3.8,2.5,2.7,3,2.6,2.2,2.1,4.2,5.8,6.2,6.3,6.5,6.6,6.3,6.2,6,5.7]                  # Load_VPP1

P_Wmax_2=[4.7,5.1,4.3,4.1,3.8,3.9,4,5,5,4.8,3.9,4.3,5,5.2,5.8,5.6,1.6,0.9,5.8,4.1,3.6,3.5,3.1,3.8]                # Wind_VPP2
load_2=[5,4,4,4.2,4.1,3.6,3.4,3.7,3.9,3.8,3.9,4,4.1,4.2,3.7,3,5.1,6.1,5.8,6.2,6.3,5.5,5,3.8]                      # Load_VPP2

P_Wmax_3=[9.3,10.1,7.2,7.5,7.9,6.4,7.1,6.9,5.6,5.4,5.2,4,3.8,3,2.8,3.2,2.5,1.1,2.1,2.9,2.7,3,4.6,5.5]             # Wind_VPP3
load_3=[4,2.1,1.1,1.1,0.7,1,1.9,3.6,3.8,4.2,5.8,5.6,5.8,5.6,5.7,6.1,8,10,9.4,8.2,6.2,5.5,4.8,2.2]                 # Load_VPP3


Load_fix=P_Wmax_1+P_Wmax_2+P_Wmax_3               # total wind power of three VPPs
Total_load=load_1+load_2+load_3                     # total load of three VPP
3*Total_load-Load_fix

λ_o= [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.75,0.75,0.75,0.75,1.2,1.2,1.2,0.75,0.75,0.75,0.75,1.2,1.2,1.2,1.2,0.4,0.4]  # energy prices, coming from the market outcome without TOU

T=24                      # length of market horizon


#-------------------------definations of lower-level problem: System operation-------------
@variable(Lower(model), 0<=P_g₁[1:T]<=P_G_max[1])                               # bounds for output of g₁
@variable(Lower(model), 0<=P_g₂[1:T]<=P_G_max[2])                               # bounds for output of g₂
@variable(Lower(model), 0<=P_g₃[1:T]<=P_G_max[3])                               # bounds for output of g₃




@constraint(Lower(model), P_g₁ +P_g₂ +P_g₃ == Total_load )                     # power balance constraint

@objective(Lower(model), Min, sum(c_G[1].*P_g₁) +sum(c_G[2].*P_g₂) +sum(c_G[3].*P_g₃)  )  # obj of lower-level problem: min generation cost


# definations and constraints of lupper-level problem: max profits of DSO   (big-M representing the indicator judgement)
#      let's double check the indicator judgement, next we use big-M to rewrite it:
#      first, the upper-level obj has this judgement(here the time interval is only 1 hour, we dropped off the whole time domain: 24 hours):
#               P_DSO=sum(P_j_VPPp-P_j_VPPs)

#             if P_DSO>=0,   which means that DSO should buy energy from PM to support VPPs operation
#                then P_DSO,p=P_DSO
#             else if P_DSO<0                               （judgement 1）
#                then P_DSO,p=0

#             if P_DSO<0,   which means that DSO should sell energy to PM to support VPPs operation
#                then P_DSO,p=-P_DSO
#             else if P_DSO>0                                (judgement 2）
#                then P_DSO,s=0
@variable(Upper(model), P_DSO_p[1:T])                                     # P_DSO_p
@variable(Upper(model), P_DSO_s[1:T])                                     # P_DSO_s
@variable(Upper(model), λ_PMs[i]<=λ_VPPp[i in 1:T]<=λ_PMp[i])             # bounds for VPP buying energy from DSO
@variable(Upper(model), λ_PMs[i]<=λ_VPPs[i in 1:T]<=λ_PMp[i])             # bounds for VPP selling energy to DSO

@variable(Upper(model), N_1[1:T] ,Bin)                                    # judgement 1
@variable(Upper(model), N_2[1:T] ,Bin)                                    # judgement 2
@constraint(Upper(model),N_1-N_2==0)                                      # N_2 & N_1 must be the same

# big-M formulation of judgement 2
@constraint(Upper(model),P_VPP_1+P_VPP_2+P_VPP_3<=M_2*N_2)
@constraint(Upper(model),-M_2*(ones(1,T)'-N_2).<=P_VPP_1+P_VPP_2+P_VPP_3)
@constraint(Upper(model),-M_2*N_2.<=P_DSO_s+(P_VPP_1+P_VPP_2+P_VPP_3))
@constraint(Upper(model),P_DSO_s+(P_VPP_1+P_VPP_2+P_VPP_3).<=M_2*N_2)
@constraint(Upper(model),-M_2*(ones(1,T)'-N_2).<=P_DSO_s)
@constraint(Upper(model),P_DSO_s.<=M_2*(ones(1,T)'-N_2))

# big-M formulation of judgement 1
@constraint(Upper(model),P_VPP_1+P_VPP_2+P_VPP_3<=M_1*N_1)
@constraint(Upper(model),-M_1*(ones(1,T)'-N_1).<=P_VPP_1+P_VPP_2+P_VPP_3)
@constraint(Upper(model),-M_1*(ones(1,T)'-N_1).<=P_DSO_p-(P_VPP_1+P_VPP_2+P_VPP_3))
@constraint(Upper(model),P_DSO_p-(P_VPP_1+P_VPP_2+P_VPP_3).<=M_1*(ones(1,T)'-N_1))
@constraint(Upper(model),-M_1*N_1.<=P_DSO_p)
@constraint(Upper(model),P_DSO_p.<=M_1*N_1)



# objs of upper-level problems and lower-level problems
@objective(Upper(model), Max, sum(λ_PMs.*P_DSO_s-λ_PMp.*P_DSO_p) + sum(λ_VPPp.*(P_VPP_1_p+P_VPP_2_p+P_VPP_3_p))-sum(λ_VPPs.*(P_VPP_1_s+P_VPP_2_s+P_VPP_3_s)))

C_VPP1=sum(λ_VPPp.*P_VPP_1_p-λ_VPPs.*P_VPP_1_s)+
a_MT[1]*sum(P_MT_1.*P_MT_1)+b_MT[1]*sum(P_MT_1)+c_MT[1]+
λ_BS[1]*sum(P_BS_1.*P_BS_1)

C_VPP2=sum(λ_VPPp.*P_VPP_2_p-λ_VPPs.*P_VPP_2_s)+
a_MT[2]*sum(P_MT_2.*P_MT_2)+b_MT[2]*sum(P_MT_2)+c_MT[2]+
λ_BS[2]*sum(P_BS_2.*P_BS_2)

C_VPP3=sum(λ_VPPp.*P_VPP_3_p-λ_VPPs.*P_VPP_3_s)+
a_MT[3]*sum(P_MT_3.*P_MT_3)+b_MT[3]*sum(P_MT_3)+c_MT[3]+
λ_BS[3]*sum(P_BS_3.*P_BS_3)

@objective(Lower(model), Min, C_VPP1+C_VPP2+C_VPP3)

#---------------------------solve the bi-level model-------------------

BilevelJuMP. set_mode(model , BilevelJuMP. StrongDualityMode())
# set_optimizer(model , SCIP.Optimizer)
# set_attribute(model, "limits/gap", 0.0280)
# set_time_limit_sec(model, 700.0)
optimize!(model)