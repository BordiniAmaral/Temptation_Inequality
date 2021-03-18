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
gridq = np.array([0, 1000, 2000, 3000, 4000, 5000])*12.0

# Starting a Identity weights matrix
W = np.identity(len(gridq))

# Building which parameters will be explored:
grid_x0 = np.arange(170,250,2)*12
grid_gam1 = np.arange(-0.2,0.7,0.02)
grid_gam2 = np.arange(0.65,0.8,0.005)

# Running GMM with all data
sqe = cb.run_x0_simulations(grid_gam1, grid_gam2, grid_x0, data_x, data_y, gridq, W, quantile)
sol = np.argwhere(sqe == np.min(sqe))
print("\nParameters found: \n x0:",grid_x0[sol[0][0]]/12,"\n gam1:",grid_gam1[sol[0][1]],"\n gam2:",grid_gam2[sol[0][2]],"\n---------------------------")

# Running bootsrap:
samples = 300
sample_size = 58034
verbose = False

boot_estimates, corner_solution, IC95_x0, IC95_gam1, IC95_gam2 = cb.bootstrap(samples, sample_size, grid_x0, grid_gam1, grid_gam2, data_x, data_y, gridq, W, quantile, verbose)

#%% Simulating

import OLG_A_CRRA as olg

# Example gridspace...
grida_space = np.zeros(600)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 50
    elif a < 12:
        grida_space[a] = grida_space[a-1]+50
    elif a < 106:
        grida_space[a] = grida_space[a-1]+100
    else:
        grida_space[a] = grida_space[a-1]*1.01

grida = np.concatenate((-np.flip(grida_space[0:(np.int(len(grida_space)/150))]),[0],grida_space))

gridz = values*1
Pi = transition

gridz[gridz<2] = 0.01 # Setting it to a close-to-zero value because zero consumption is reserved to identify unfeasible scenarios.

n = 40 # always verify if size is adequate for gridz size
xi = 1.072
delta = 0.035
beta = 0.96
alpha = 0.33
sigma_x = 0.21
sigma_y = 0.29
transfer = 0
x0 = 210*12

A = 1 
r_low = 0.04
r_high = 0.1
# Check if using correct mass grid for z (it must be coherent with gridz imported from PNAD)
mass_z = np.concatenate((np.repeat(0.1,9),np.repeat(0.02,5))) 

# print("Elapsed (with compilation) = %s" % (end - start))# 1.Computing baseline GE: No temptation, no function T
start = time.time()
KL1, w1, C1, x1, y1, V1, choice_a1, distr_mass1, k_mass1, k_total1, c_mass1, c_total1, r1, init_r1 = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi = 0, r_low = r_low, r_high = r_high, mass_z = mass_z, transfer = transfer, A = A, x0 = 0, temptation = True, tol = 1e-2, maxit = 10000)
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))

# 2.Computing GE with non-homothetic preferences, matching previous aggregate savings
step = 0.01
tol = 1e-3
temptation = False
beta_equivalent2, results_equivalent2 = olg.ge_match_by_capital(beta, step, tol, k_total1, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r1, mass_z, transfer, A, x0, temptation)

choice_a_eq2 = results_equivalent2[6]
distr_mass_eq2 = results_equivalent2[7]
k_mass_eq2 = results_equivalent2[8]

# 3.Computing GE with full temptation preferences, matching previous aggregate savings
step = 0.01
tol = 1e-3
temptation = True
beta_equivalent3, results_equivalent3 = olg.ge_match_by_capital(beta, step, tol, k_total1, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r1, mass_z, transfer, A, x0, temptation)

choice_a_eq3 = results_equivalent3[6]
distr_mass_eq3 = results_equivalent3[7]
k_mass_eq3 = results_equivalent3[8]

# Results equivalent indexes:
    #  0  1  2  3  4  5      6         7         8       9       10       11     12
    # KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, r

#%% Plotting

import StatsAndGraphs as stg

# Plotting comparisons
mass_by_k1, mass_by_age_k1 = stg.capital_distributions(n, grida, gridz, distr_mass1, k_mass1)
mass_by_k_eq2, mass_by_age_k_eq2 = stg.capital_distributions(n, grida, gridz, distr_mass_eq2, k_mass_eq2)
mass_by_k_eq3, mass_by_age_k_eq3 = stg.capital_distributions(n, grida, gridz, distr_mass_eq3, k_mass_eq3)

stg.compare_total_k_distr(mass_by_k1, mass_by_k_eq2, mass_by_k_eq3, grida, bin_size = 2500, description = "Calibration (B)", label1 = "Baseline CRRA", label2 = "Non-homothetic", label3 = "Temptation", log = False, trim_upper = True, trim_value = 600000)
stg.compare_total_k_lorenz(mass_by_k1, mass_by_k_eq2, mass_by_k_eq3, grida, description = "Calibration (B)", label1 = "Baseline CRRA", label2 = "Non-homothetic", label3 = "Temptation")

# Plotting Asset evolution by quantile
age_start = 25
quants = np.concatenate((np.repeat(0.1,9),np.repeat(0.02,5))) 

stg.plot_k_evolution(age_start, n, mass_by_age_k1, quants, grida, description = "Baseline CRRA")
stg.plot_k_evolution(age_start, n, mass_by_age_k_eq2, quants, grida, description = "Non homothetic")
stg.plot_k_evolution(age_start, n, mass_by_age_k_eq3, quants, grida, description = "Temptation")

# Savings Rate

# Default quants
quants = np.concatenate((np.arange(0,1,0.1),np.arange(0.92,1.01,0.02))) 
# Alternative
# quants = np.arange(0,1.01,0.25)

# Checking some stats
include_interest = True
quant_mean_eq, quant_wt_sd_eq = stg.savings_by_quants(n, grida, choice_a_eq, gridz, r_eq, distr_mass_eq, quants, include_interest)
quant_mean_nt, quant_wt_sd_nt = stg.savings_by_quants(n, grida, choice_a_nt, gridz, r_nt, distr_mass_nt, quants, include_interest)

# Comparing lifecycle average savings
age_start = 25
quants = np.array([0,0.5,0.9,0.99,1])
include_interest = True
n_select = n-5 # Dissavings at the end of life are huge, making visualization poor

quant_value_nt, quant_value_eq = stg.compare_savings_rate(age_start, n, n_select, quants, grida, gridz, r_nt, r_eq, choice_a_nt, choice_a_eq, distr_mass_nt, distr_mass_eq, include_interest, description1 = "No Temptation", description2 = "With Temptation")
quants = np.array([0.5,0.9,0.99])
stg.compare_k_evolution(age_start,n,mass_by_age_k_nt, mass_by_age_k_eq, quants, grida, description1 = "No Temptation", description2 = "With Temptation")

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
olg.calculate_wealth_gini(n, grida, gridz, distr_mass_eq2, k_mass_eq2)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass_eq3, k_mass_eq3)

# Wealth report
show_zero = True
quants = np.array([0.5,0.9, 0.95, 0.99])

stg.savings_and_wealth_report(n, mass_by_k1, grida, quants, distr_mass1, show_zero)
stg.savings_and_wealth_report(n, mass_by_k_eq2, grida, quants, distr_mass_eq2, show_zero)
stg.savings_and_wealth_report(n, mass_by_k_eq3, grida, quants, distr_mass_eq3, show_zero)