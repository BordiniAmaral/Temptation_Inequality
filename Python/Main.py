##############################################################################
# Project:      Masters Dissertation - Quantitative Heterogeneous Temptation #
# Author:       Felipe Bordini do Amaral                                     #
# Institution:  FGV - EESP                                                   #
# Year:         2020                                                         #
#                                                                            #
# Oriented by:  Prof. Tiago Cavalcanti                                       #
#               Prof. Pierluca Pannella                                      #
#                                                                            #
##############################################################################
#                               MAIN FILE                                    #
##############################################################################

#%% Loading some packages

import pandas as pd
import numpy as np
import time

#%% Income Process

import Income_Process

import_path = "D:\Google Drive\Working Cloud\EESP-Mestrado\Dissertação\PNADC\R Project\Transitions by age"
files = ["average_matrix1.csv", "quantile_frontiers.csv", "quantile_values.csv"]

transition, values = Income_Process.load_income_process(import_path, files)

#%% Calibration

import Calibration as cb

## ONLY POF ANALYSIS

# Loading empirical data

pof_path = "D:\Google Drive\Working Cloud\EESP-Mestrado\Dissertação\POF\R Project\Exports\\temptation_filtered.csv"
pof_df = cb.import_POF_data(pof_path)
data_x, data_y = cb.select_data(pof_df)

# Selecting moment grid (if not quantile, using annual consumption per capita grid)
# Remember: include zero, leave upper grid open

# quantile = True # Pick this if selecting moments by sample quantile
# gridq = np.array([0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99])

quantile = False # Pick this if selecting moments by arbitrary total consumption levels
gridq = np.array([0, 500, 1000, 2000, 3000, 4000, 5000])*12.0

# Starting a Identity weights matrix
W = np.identity(len(gridq))

# Building which parameters will be explored:
grid_x0 = np.arange(110,160,1)*12
grid_gam1 = np.arange(-1, 0, 0.02)
grid_gam2 = np.arange(0.72,0.85,0.005)

# Running GMM with all data
sqe = cb.run_x0_simulations(grid_gam1, grid_gam2, grid_x0, data_x, data_y, gridq, W, quantile,verbose = True)
sol = np.argwhere(sqe == np.min(sqe))
print("\nParameters found: \n x0:",grid_x0[sol[0][0]]/12,"\n gam1:",grid_gam1[sol[0][1]],"\n gam2:",grid_gam2[sol[0][2]],"\n---------------------------")

# Running bootsrap:
samples = 50
sample_size = 58034
verbose = False

boot_estimates, corner_solution, IC95_x0, IC95_gam1, IC95_gam2 = cb.bootstrap(samples, sample_size, grid_x0, grid_gam1, grid_gam2, data_x, data_y, gridq, W, quantile, verbose)

# In case you want to save bootstrap results
boot_df = pd.DataFrame(boot_estimates, index = None, columns = ["x0","gam1","gam2"])
sim_id = str("003_02")
boot_df.to_csv("D:\Google Drive\Working Cloud\EESP-Mestrado\Dissertação\Temptation\Exports\Python (Post Jan21)\sim_"+sim_id+"_bootstrap.csv", sep = ";")

# In case you want to load bootstrap results
sim_id = str("003_01")
boot_df = pd.read_csv("D:\Google Drive\Working Cloud\EESP-Mestrado\Dissertação\Temptation\Exports\Python (Post Jan21)\sim_"+sim_id+"_bootstrap.csv", sep = ";")

#%% Simulating

import OLG_A_CRRA as olg

# Example gridspace...
grida_space = np.zeros(750)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 30
    elif a < 8:
        grida_space[a] = grida_space[a-1]+30
    else:
        grida_space[a] = grida_space[a-1]*1.011

grida = np.concatenate((-np.flip(grida_space[0:(np.int(len(grida_space)/150))]),[0],grida_space))

'''
# Smaller (and quicker) gridspace for testing

grida_space = np.zeros(220)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 50
    elif a == 1:
        grida_space[a] = 100
    elif a < 8:
        grida_space[a] = grida_space[a-1]*1.4
    elif a < 12:
        grida_space[a] = grida_space[a-1]*1.3
    elif a < 15:
        grida_space[a] = grida_space[a-1]*1.18
    elif a < 20:
        grida_space[a] = grida_space[a-1]*1.09
    else:
        grida_space[a] = grida_space[a-1]*1.025

grida = np.concatenate((-np.flip(grida_space[0:(np.int(len(grida_space)/20))]),[0],grida_space))


# Grid with a "high" minimum investment

grida_space = np.zeros(460)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 1200
    elif a < 106:
        grida_space[a] = grida_space[a-1]*1.005
    elif a < 200:
        grida_space[a] = grida_space[a-1]*1.01
    else:
        grida_space[a] = grida_space[a-1]*1.02

grida = np.concatenate(([0],grida_space))

'''

# Check if using correct mass grid for z (it must be coherent with gridz imported from PNAD)
mass_z = np.concatenate((np.repeat(0.1,9),np.repeat(0.02,5))) 
gridz = values*1

# Using PNAD calibration:
Pi = transition*1

# Using arbitrary transition:
# Pi = Income_Process.load_test_transition(mass_z, persist = 1, width = 1)

gridz[gridz<2] = 0.01 # Setting it to a close-to-zero value because zero consumption is reserved to identify unfeasible scenarios.

n = 40 # always verify if size is adequate for gridz size
xi =  0.2494 # 0.2494
delta = 0.035
beta = 0.96
alpha = 0.33
sigma_x = 2.17 # 2.17
sigma_y = 2.67 # 2.67 
transfer = 0
x0 = 133*12
x0_1 = 0
xi_1 = 0.2494 # 0.2494
inv_cost = 0

A = 1 
r_low = 0.04
r_high = 0.1

#---------------------------------------- Benchmark CRRA -------------------------------------------

# 1.Computing baseline GE: No temptation, no minimum consumption
start_time = time.time()
KL1, w1, C1, x1, y1, V1, choice_a1, distr_mass1, k_mass1, k_total1, c_opt1, r1, init_r1, U1, T1 = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi = xi_1, r_low = r_low, r_high = r_high, mass_z = mass_z, transfer = transfer, A = A, x0 = x0_1, inv_cost = inv_cost, temptation = False, tol = 1e-2, maxit = 10000)
end_time = time.time()
hour = round(((end_time - start_time)/60)//60)
minute = round((end_time - start_time)//60 - hour*60)
second = round((end_time - start_time) - hour*60*60 - minute*60)
print("Finished! \nTime elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")

#-------------------------------- Exercises with recalibrated Beta ---------------------------------

# 2a.Computing GE with minimum consumption, but not temptation, matching previous aggregate savings
step = 0.01
tol = 5e-3
temptation = False
beta_equivalent2, results_equivalent2 = olg.ge_match_by_capital(beta, step, tol, k_total1, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi_1, r1, mass_z, transfer, A, x0, temptation, inv_cost)

V_eq2 = results_equivalent2[5]
choice_a_eq2 = results_equivalent2[6]
distr_mass_eq2 = results_equivalent2[7]
k_mass_eq2 = results_equivalent2[8]
k_total_eq2 = results_equivalent2[9]
c_opt_eq2 = results_equivalent2[10]

# 3a. Computing GE with temptation, but not minimum consumption, matching previous aggregate savings
step = 0.01
tol = 5e-3
temptation = True
beta_equivalent3, results_equivalent3 = olg.ge_match_by_capital(beta, step, tol, k_total1, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r1, mass_z, transfer, A, x0_1, temptation, inv_cost)

V_eq3 = results_equivalent3[5]
choice_a_eq3 = results_equivalent3[6]
distr_mass_eq3 = results_equivalent3[7]
k_mass_eq3 = results_equivalent3[8]
k_total_eq3 = results_equivalent3[9]
c_opt_eq3 = results_equivalent3[10]

'''
# 4a.Computing GE with temptation and minimum consumption, matching previous aggregate savings
step = 0.01
tol = 5e-3
temptation = True
beta_equivalent4, results_equivalent4 = olg.ge_match_by_capital(beta, step, tol, k_total1, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r1, mass_z, transfer, A, x0, temptation, inv_cost)

V_eq4 = results_equivalent4[5]
choice_a_eq4 = results_equivalent4[6]
distr_mass_eq4 = results_equivalent4[7]
k_mass_eq4 = results_equivalent4[8]
k_total_eq4 = results_equivalent4[9]
c_opt_eq3 = results_equivalent3[10]
'''

# Results equivalent indexes:
    #  0  1  2  3  4  5      6         7         8       9       10    11  
    # KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_opt,  r
    
#-------------------------------- Exercises with constant Beta ---------------------------------

# 2b.Computing GE with minimum consumption, but not temptation, matching previous aggregate savings
temptation = False
start_time = time.time()
KL2, w2, C2, x2, y2, V2, choice_a2, distr_mass2, k_mass2, k_total2, c_opt2, r2, init_r2, U2, T2 = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi = xi_1, r_low = r_low, r_high = r_high, mass_z = mass_z, transfer = transfer, A = A, x0 = x0, inv_cost = inv_cost, temptation = temptation, tol = 1e-2, maxit = 10000)
end_time = time.time()
hour = round(((end_time - start_time)/60)//60)
minute = round((end_time - start_time)//60 - hour*60)
second = round((end_time - start_time) - hour*60*60 - minute*60)
print("Finished! \nTime elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")

# 3b.Computing GE with temptation, but not minimum consumption
temptation = True
start_time = time.time()
KL3, w3, C3, x3, y3, V3, choice_a3, distr_mass3, k_mass3, k_total3, c_opt3, r3, init_r3, U3, T3 = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi = xi, r_low = r_low, r_high = r_high, mass_z = mass_z, transfer = transfer, A = A, x0 = x0_1, inv_cost = inv_cost, temptation = temptation, tol = 1e-2, maxit = 10000)
end_time = time.time()
hour = round(((end_time - start_time)/60)//60)
minute = round((end_time - start_time)//60 - hour*60)
second = round((end_time - start_time) - hour*60*60 - minute*60)
print("Finished! \nTime elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")

'''
# 4b.Computing GE with temptation and minimum consumption
temptation = True
start_time = time.time()
KL4, w4, C4, x4, y4, V4, choice_a4, distr_mass4, k_mass4, k_total4, c_opt4, r4, init_r4, U4, T4 = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0, inv_cost, temptation, tol = 1e-2, maxit = 10000)
end_time = time.time()
hour = round(((end_time - start_time)/60)//60)
minute = round((end_time - start_time)//60 - hour*60)
second = round((end_time - start_time) - hour*60*60 - minute*60)
print("Finished! \nTime elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")
'''
#%% Plotting

import StatsAndGraphs as stg

# Plotting comparisons
mass_by_k1, mass_by_age_k1 = stg.capital_distributions(n, grida, gridz, distr_mass1, k_mass1)
mass_by_k2, mass_by_age_k2 = stg.capital_distributions(n, grida, gridz, distr_mass2, k_mass2)
mass_by_k3, mass_by_age_k3 = stg.capital_distributions(n, grida, gridz, distr_mass3, k_mass3)
#mass_by_k4, mass_by_age_k4 = stg.capital_distributions(n, grida, gridz, distr_mass4, k_mass4)
mass_by_k_eq2, mass_by_age_k_eq2 = stg.capital_distributions(n, grida, gridz, distr_mass_eq2, k_mass_eq2)
mass_by_k_eq3, mass_by_age_k_eq3 = stg.capital_distributions(n, grida, gridz, distr_mass_eq3, k_mass_eq3)
#mass_by_k_eq4, mass_by_age_k_eq4 = stg.capital_distributions(n, grida, gridz, distr_mass_eq4, k_mass_eq4)

# zero_plot = np.zeros(len(grida))
# 1 - CRRA only
# 2 - Stone-Geary
# 3 - Temptation
# 4 - Stone-Geary + Temptation
trim_value = 60000
description = "Fixed Aggregate Savings"
mass1 = mass_by_k1
label1 = "Baseline CRRA"
mass2 = mass_by_k_eq2
label2 = "Minimum Consumption"
mass3 = mass_by_k_eq3
label3 = "Temptation"
stg.compare_total_k_distr(mass1, mass2, mass3, grida, bin_size = 3000, description = description, label1 = label1, label2 = label2, label3 = label3, log = False, trim_upper = True, trim_value = trim_value)
stg.compare_total_k_lorenz(mass1, mass2, mass3, grida, description, label1, label2, label3)

# Comparing lifecycle average savings
age_start = 25
quants = np.array([0,0.5,0.9,0.99,1])
include_interest = True
n_select = n-5 # Dissavings at the end of life are huge, making visualization poor

quant_value_nt, quant_value_eq = stg.compare_savings_rate(age_start, n, n_select, quants, grida, gridz, \
                                                          r1, r2, r3, \
                                                          choice_a1, choice_a2, choice_a3, \
                                                          distr_mass1, distr_mass2, distr_mass3, \
                                                          include_interest, \
                                                          description1 = "Baseline CRRA", description2 = "Minimum Consumption", description3 = "Temptation")


quants = np.array([0.5,0.9,0.99])
stg.compare_k_evolution(age_start,n,mass_by_age_k1, mass_by_age_k_eq2, mass_by_age_k_eq3, quants, grida, description1 = "Baseline CRRA", description2 = "Non-homothetic", description3 = "Temptation")

# Exporting savings rate curves
savings_nt = pd.DataFrame(quant_value_nt, index = quants, columns = np.arange(age_start,age_start+n))
savings_eq = pd.DataFrame(quant_value_eq, index = quants, columns = np.arange(age_start,age_start+n))

# Insert simulation ID to find this result later
sim_id = str("xxx")
savings_nt.to_csv("D:\Google Drive\Working Cloud\EESP-Mestrado\Dissertação\Temptation\Exports\Python (Post Jan21)\sim_"+sim_id+"_savings_nt.csv", sep = ";")
savings_eq.to_csv("D:\Google Drive\Working Cloud\EESP-Mestrado\Dissertação\Temptation\Exports\Python (Post Jan21)\sim_"+sim_id+"_savings_eq.csv", sep = ";")

#%% Calculating some stats

# Wealth Gini
olg.calculate_wealth_gini(n, grida, gridz, distr_mass1, k_mass1)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass2, k_mass2)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass3, k_mass3)
#olg.calculate_wealth_gini(n, grida, gridz, distr_mass4, k_mass4)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass_eq2, k_mass_eq2)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass_eq3, k_mass_eq3)
#olg.calculate_wealth_gini(n, grida, gridz, distr_mass_eq4, k_mass_eq4)

# Wealth report
show_zero = True
quants = np.array([0.5,0.9, 0.95, 0.99])

stg.savings_and_wealth_report(n, mass_by_k1, grida, quants, distr_mass1, show_zero)
stg.savings_and_wealth_report(n, mass_by_k2, grida, quants, distr_mass2, show_zero)
stg.savings_and_wealth_report(n, mass_by_k3, grida, quants, distr_mass3, show_zero)
# stg.savings_and_wealth_report(n, mass_by_k4, grida, quants, distr_mass4, show_zero)
stg.savings_and_wealth_report(n, mass_by_k_eq2, grida, quants, distr_mass_eq2, show_zero)
stg.savings_and_wealth_report(n, mass_by_k_eq3, grida, quants, distr_mass_eq3, show_zero)
# stg.savings_and_wealth_report(n, mass_by_k_eq4, grida, quants, distr_mass_eq4, show_zero)

# Simulated fractions
gam1 = -0.52
gam2 = 0.794
verbose = True
quantile = False # Pick this if selecting moments by arbitrary total consumption levels
gridq = np.array([0, 500, 1000, 2000, 3000, 4000, 5000])*12.0
fraction_sim = cb.compute_sim_tempt_frac(data_x, data_y, x0, gam1, gam2, gridq, quantile, verbose)

# Capital supply and demand curves

grid_r = np.arange(0.037,0.055,0.001)
k_supply1 = olg.compute_k_supply_curve(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi_1, grid_r, r1, mass_z, transfer, A, temptation = False, x0 = x0_1)
k_supply2 = olg.compute_k_supply_curve(n, beta_equivalent2, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi_1, grid_r, r1, mass_z, transfer, A, temptation = False, x0 = x0)
k_supply3 = olg.compute_k_supply_curve(n, beta_equivalent3, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, grid_r, r1, mass_z, transfer, A, temptation = True, x0 = x0_1)

k_demand = olg.compute_k_demand_curve(grid_r, delta, w1, A, gridz, mass_z, n, alpha)

stg.plot_k_curves(k_supply1 ,k_supply2, k_supply3, grid_r, k_demand, label1 = "Baseline CRRA", label2 = "Minimum Consumption", label3 = "Temptation")

# ------------------------------ Rascunhos -----------------------------------

# Getting dxdc graph to determine temptation tax

import matplotlib.pyplot as plt
import Temptation_CRRA as tpt

# Adjust the parameters you want to explore here
xi1 = 0.2494
sigma_x1 = 2.17
sigma_y1 = 2.67
x01 = 133*12

xi2 = 2.7e5
sigma_x2 = 2.17
sigma_y2 = 4.34
x02 = 133*12

C = np.arange(134*12,3000*12,1*12)*1.0
x1, y1, U1, T1, frac1, dxdc1, dudc1, dtdc1, dudx1, dtdy1 = tpt.calculate_allocation_by_C(C, x01, sigma_x1, sigma_y1, xi1)
x2, y2, U2, T2, frac2, dxdc2, dudc2, dtdc2, dudx2, dtdy2 = tpt.calculate_allocation_by_C(C, x02, sigma_x2, sigma_y2, xi2)

plt.figure(figsize=(7,7))
plt.plot(C/12,dxdc1, label = 'dx(c)/dc (calibration)')
plt.plot(C/12,dxdc2, label = 'dx(c)/dc (alternative)')
plt.ylabel('dx(c)/dc')
plt.xlabel('C.onsumption (R$ montlhy p.c.)')
plt.legend()
plt.ylim([0.6, 1])
plt.title('Temptation intertemporal tax')


plt.figure(figsize=(7,7))
plt.plot(C[80:120]/12, (dudc+dtdc)[80:120], label = "Period t      : [d(U+T)/dc] (c)")
plt.plot(C[80:120]/12, beta*dudc[80:120], label = "Period t+1 : \u03B2[dU(x)/dx * dx(c)/dc] (c)", linestyle='--')
plt.legend()
plt.title('Consumption today vs. Consumption tomorrow')
plt.xlabel('Total Consumption (R$ montlhy p.c.)')
