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




# ==============================================
# 24-hour TOU Load-Price Elasticity Matrix Model
# Author: (Your Name)
# ==============================================

using LinearAlgebra

# --- Define 24 Time Periods ---
n_periods = 24

# Define TOU period type for each hour: 1=valley, 2=flat, 3=peak
tou_type = [
    1,1,1,1,1,  # 0–5h: valley
    2,2,2,      # 6–8h: flat
    3,3,3,3,    # 9–12h: peak
    2,2,2,2,    # 13–16h: flat
    3,3,3,      # 17–19h: peak
    2,2,2,      # 20–22h: flat
    1           # 23h: valley
]

# --- Base Load (kWh) for each hour ---
base_load = [
    800, 750, 700, 680, 700,   # 0–5h valley
    1000, 1200, 1400,          # 6–8h flat
    2000, 2200, 2300, 2400,    # 9–12h peak
    1800, 1600, 1500, 1400,    # 13–16h flat
    2500, 2600, 2400,          # 17–19h peak
    1700, 1500, 1300,          # 20–22h flat
    900                        # 23h valley
]

# --- Base TOU Prices (¥/kWh) ---
price_valley = 0.3
price_flat   = 0.6
price_peak   = 1.0

base_price = [price_valley, price_flat, price_peak][tou_type]

# --- Define Elasticity Parameters ---
self_elasticity = -0.25       # own-price elasticity
cross_elasticity_adj = 0.05   # adjacent-period elasticity
cross_elasticity_nonadj = 0.02 # non-adjacent weak elasticity

# --- Construct Elasticity Matrix (24x24) ---
E = zeros(Float64, n_periods, n_periods)
for i in 1:n_periods
    for j in 1:n_periods
        if i == j
            E[i, j] = self_elasticity
        elseif abs(i - j) == 1
            E[i, j] = cross_elasticity_adj
        else
            E[i, j] = cross_elasticity_nonadj
        end
    end
end

# --- Define New TOU Prices (¥/kWh) ---
# Example: slightly increase peak prices, reduce valley prices
new_price = copy(base_price)
for i in 1:n_periods
    if tou_type[i] == 3      # peak
        new_price[i] *= 1.10
    elseif tou_type[i] == 1  # valley
        new_price[i] *= 0.90
    end
end

# --- Compute Load Change ---
ΔP = (new_price .- base_price) ./ base_price
ΔL = E * ΔP
new_load = base_load .* (1 .+ ΔL)

# --- Display Results ---
println("=== TOU Elasticity Matrix Model ===")
println("Total Base Load (kWh): ", sum(base_load))
println("Total New Load (kWh): ", round(sum(new_load), digits=2))
println("Total Load Change: ", round(sum(new_load) - sum(base_load), digits=2), " kWh\n")

println("Hour | TOU | Base_Price | New_Price | Base_Load | New_Load | ΔLoad(%)")
for i in 1:n_periods
    println(rpad(i,4), " | ",
            tou_type[i], "   | ",
            round(base_price[i], digits=2), "       | ",
            round(new_price[i], digits=2), "     | ",
            round(base_load[i], digits=0), "      | ",
            round(new_load[i], digits=0), "    | ",
            round(ΔL[i]*100, digits=2))
end

# ==============================================
# Notes:
# - Negative self-elasticity: load drops when price rises.
# - Positive cross-elasticity: load shifts between adjacent hours.
# - You can modify elasticity parameters for sensitivity analysis.
# ==============================================



%function fxx=DR(fx) %负荷响应
%% 调这里的参数就行
fx=p1';%负荷曲线 (单独初始负荷或者汽车   or 初始负荷+汽车)
fj=1.18; %峰时段电价
pj=0.72; %平时段电价
gj=0.43; %谷时段电价
DIANJIA=[gj,gj,gj,gj,pj,pj,pj,pj,pj,pj,fj,fj,fj,fj,pj,pj,pj,fj,fj,fj,fj,fj,gj,gj];% 各时刻电价 24h

El=-0.3;%自弹性系数
Efp=0.03;%峰-平弹性系数
Efg=0.05;%峰-谷弹性系数
Epg=0.03;%平-谷弹性系数

%%
kf=(fj-pj)/pj;
kp=0;
kg=(gj-pj)/pj;

lam=zeros(24,24);
for i=1:24
    if DIANJIA(i)==fj
        lam(i,i)=kf*El;
    elseif DIANJIA(i)==pj
        lam(i,i)=0;
    elseif DIANJIA(i)==gj
        lam(i,i)=kg*El;
    end
    for j=i+1:24
        if DIANJIA(i)==fj && DIANJIA(j)==gj
        lam(i,j)=kg*Efg;
        elseif DIANJIA(i)==fj && DIANJIA(j)==pj
            lam(i,j)=-(kf-kg)*Efp;
        elseif DIANJIA(i)==pj && DIANJIA(j)==gj
            lam(i,j)=kg*Efp;
        elseif DIANJIA(i)==gj && DIANJIA(j)==fj
        lam(i,j)=-kg*Efg;
        elseif DIANJIA(i)==pj && DIANJIA(j)==gj
            lam(i,j)=(kf-kg)*Efp;
        elseif DIANJIA(i)==gj && DIANJIA(j)==pj
            lam(i,j)=-kg*Efp;
        end
        lam(j,i)=-lam(i,j);
    end
end
fxx=fx+lam*fx;
figure;
plot(fxx+p0','r-o')
hold on
plot(p0'+p1','b-*')
legend('需求响应后','需求响应前');
xlabel('时间');
ylabel('功率');
axis([1 24 0 3500]);
hold on
yyaxis right
%x=[0.5  2.2 3.2 4.2 5.2 6.2 7.2 8.2 9.2 10.2 11.2 12.2 13.2 14.2 15.2 16.2 17.2 18.2 19.2 20.2 21.2 22.2 23.2 24.5];
%x=linspace(0.5, 24, 24);
% x=[0 1 2 3 4 5.5 7 8 9 10 11 12 13 14 15 16 16.5 18 19.4 20 21 21.5 23.2 25];
x=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
stairs(x,DIANJIA,'b-')
xlabel('时间(h)','Fontname','宋体','Fontsize',7.5);
ylabel('电价价格(元/kWh)','Fontname','宋体','Fontsize',7.5)
axis([1 24 0 2]);
set(gca,'position',[0.09 0.13 0.8 0.65]);
set(gca,'box','off');

