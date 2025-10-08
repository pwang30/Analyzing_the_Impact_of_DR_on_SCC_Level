# Analyzing the Role of DSO in Electricity Trading of VPPs via a Stackelberg Game Model
# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Note 1: This code is full, easy to understand and replicate. Every line has its corresponding notes for a better understanding.
# Note 2: For each model of VPP/Lower level, they have same codeing structure, the only difference is the parameters. 
#         So, once you know how to structur anyone of them, please dive into DSO/Upper level modeling.

# 19.Dec.2024

#----------------------Pkgs introduction-----------------------

import Pkg 
using JuMP, BilevelJuMP,Plots,MAT, Ipopt  

#----------------------define model----------------------------

model= BilevelModel()

#----------------------relevant parameters input---------------

P_Wmax_1=[2,1.5,1.6,1.8,1.3,0.6,2.8,3.3,3.9,4,3.3,2.9,2.7,2,0.2,3.2,5.1,3.1,1.8,2,1.3,1,2,3.8]                    # Wind_VPP1
load_1=[2.2,1.8,3,6,5.8,5.2,5.6,3.8,2.5,2.7,3,2.6,2.2,2.1,4.2,5.8,6.2,6.3,6.5,6.6,6.3,6.2,6,5.7]                  # Load_VPP1

P_Wmax_2=[4.7,5.1,4.3,4.1,3.8,3.9,4,5,5,4.8,3.9,4.3,5,5.2,5.8,5.6,1.6,0.9,5.8,4.1,3.6,3.5,3.1,3.8]                # Wind_VPP2
load_2=[5,4,4,4.2,4.1,3.6,3.4,3.7,3.9,3.8,3.9,4,4.1,4.2,3.7,3,5.1,6.1,5.8,6.2,6.3,5.5,5,3.8]                      # Load_VPP2

P_Wmax_3=[9.3,10.1,7.2,7.5,7.9,6.4,7.1,6.9,5.6,5.4,5.2,4,3.8,3,2.8,3.2,2.5,1.1,2.1,2.9,2.7,3,4.6,5.5]             # Wind_VPP3
load_3=[4,2.1,1.1,1.1,0.7,1,1.9,3.6,3.8,4.2,5.8,5.6,5.8,5.6,5.7,6.1,8,10,9.4,8.2,6.2,5.5,4.8,2.2]                 # Load_VPP3


Load_fix=P_Wmax_1+P_Wmax_2+P_Wmax_3               # total wind power of three VPPs
Total_load=load_1+load_2+load_3                     # total load of three VPP

plot(3*Total_load-Load_fix)

λ_o= [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.75,0.75,0.75,0.75,1.2,1.2,1.2,0.75,0.75,0.75,0.75,1.2,1.2,1.2,1.2,0.4,0.4]  # energy prices, coming from the market outcome without TOU

T=24                      # length of market horizon


#-------------------------definations of lower-level problem: System operation-------------
@variable(Lower(model), 0<=P_g₁[1:T]<=P_G_max[1])                               # bounds for output of g₁
@variable(Lower(model), 0<=P_g₂[1:T]<=P_G_max[2])                               # bounds for output of g₂
@variable(Lower(model), 0<=P_g₃[1:T]<=P_G_max[3])                               # bounds for output of g₃

@constraint(Lower(model), P_g₁ +P_g₂ +P_g₃ == Total_load )                     # power balance constraint

@objective(Lower(model), Min, sum(c_G[1].*P_g₁) +sum(c_G[2].*P_g₂) +sum(c_G[3].*P_g₃)  )  # obj of lower-level problem: min generation cost



#-------------------------definations of upper-level problem: Strategic bidding-------------
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
