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

import numpy as np

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

# Selecting a moment grid for calibration. Options: "arbitrary" (adapt code), "percentile5", "percentile10"
# gridm = cb.create_cons_moment_grid("percentile5")

# Selecting percentile grid for moment calibration
gridq = np.arange(0,1,0.1)

# Starting a Identity weights matrix
W = np.identity(len(gridq))

# Building which parameters will be explored:
grid_x0 = np.arange(35*12,65*12,1*12)
grid_gam1 = np.arange(-2.5,-1.5,0.02)
grid_gam2 = np.arange(0.80,0.99,0.005)

# Running the simulations
sqe = cb.run_x0_simulations(grid_gam1, grid_gam2, grid_x0, data_x, data_y, gridq, W)
sol = np.argwhere(sqe == np.min(sqe))
print("\nParameters found: \n x0:",grid_x0[sol[0][0]]/12,"\n gam1:",grid_gam1[sol[0][1]],"\n gam2:",grid_gam2[sol[0][2]],"\n---------------------------")


# Running Two-Step simulations # NEED REVISION
# sqe_I, sqe_W, W_new, Omega = cb.run_two_steps(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi)

#%% Simulating

import OLG_A_CRRA as olg

# Example gridspace...
grida_space = np.zeros(200)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 100
    elif a == 1:
        grida_space[a] = 200
    elif a < 5:
        grida_space[a] = grida_space[a-1]*1.5
    elif a < 8:
        grida_space[a] = grida_space[a-1]*1.4
    elif a < 15:
        grida_space[a] = grida_space[a-1]*1.2
    elif a < 20:
        grida_space[a] = grida_space[a-1]*1.1
    else:
        grida_space[a] = grida_space[a-1]*1.03

grida = np.concatenate((-np.flip(grida_space[0:(np.int(len(grida_space)/5))]),[0],grida_space))

gridz = values
Pi = transition
n = 40
xi = 0.618
delta = 0.05
beta = 0.96
alpha = 0.4
sigma_x = 0.21
sigma_y = 0.22
transfer = 0
x0 = 45*12

A = 1 
r_low = 0.04
r_high = 0.06
# Check if using correct mass grid for z (it must be coherent with gridz imported from PNAD)
mass_z = np.concatenate((np.repeat(0.1,9),np.repeat(0.02,5))) 

# Computing a single GE - with Temptation
KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, r, init_r = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0, temptation = True, tol = 1e-2, maxit = 10000)

# Computing a single GE - without Temptation
KL_nt, w_nt, C_nt, x_nt, y_nt, V_nt, choice_a_nt, distr_mass_nt, k_mass_nt, k_total_nt, c_mass_nt, c_total_nt, r_nt, init_r_nt = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0, temptation = False, tol = 1e-2, maxit = 10000)

# Finding equivalent beta to match non-temptation aggregate capital
step = 0.01
tol = 1e-3
beta_equivalent, results_equivalent = olg.ge_match_capital(beta, step, tol, k_total_nt, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0)

#%% Plotting

import StatsAndGraphs as stg

mass_by_k, mass_by_age_k = stg.capital_distributions(n, grida, gridz, distr_mass_nt, k_mass_nt)
stg.plot_total_k_distr(mass_by_k, grida, 10000, "Calibration A, With Temptation")
stg.plot_total_k_lorenz(mass_by_k, grida, "Calibration A, With Temptation")

# Plotting comparison: Temptation vs Non-Temptation
mass_by_k, mass_by_age_k = stg.capital_distributions(n, grida, gridz, distr_mass, k_mass)
mass_by_k_nt, mass_by_age_k_nt = stg.capital_distributions(n, grida, gridz, distr_mass_nt, k_mass_nt)

stg.compare_total_k_distr(mass_by_k, mass_by_k_nt, grida, bin_size = 2500, description = "Calibration (B)", label1 = "With Temptation", label2 = "Without Temptation", log = False, trim_upper = True, trim_value = 600000)
stg.compare_total_k_lorenz(mass_by_k, mass_by_k_nt, grida, description = "Calibration (B)", label1  = "With Temptation", label2 = "Without Temptation")

# Plotting Asset evolution by quantile
age_start = 25
quants = np.concatenate((np.repeat(0.1,9),np.repeat(0.02,5))) 


stg.plot_k_evolution(age_start, n, mass_by_age_k, quants, grida, description = "with temptation")
stg.plot_k_evolution(age_start, n, mass_by_age_k_nt, quants, grida, description = "without temptation")

# Savings Rate

# Default quants
# quants = np.concatenate((np.arange(0,1,0.1),np.arange(0.92,1.01,0.02))) 
# Alternative
quants = np.arange(0,1.01,0.25)

quant_mean, quant_wt_sd = stg.savings_by_quants(n, grida, choice_a, gridz, r, distr_mass, quants)
quant_mean_nt, quant_wt_sd_nt = stg.savings_by_quants(n, grida, choice_a_nt, gridz, r_nt, distr_mass_nt, quants)

#%% Calculating some stats

# Wealth Gini
olg.calculate_wealth_gini(n, grida, gridz, distr_mass, k_mass)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass_nt, k_mass_nt)

