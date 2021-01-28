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

#%% Example Runs

# Running plain example run -------------------------------------------------
import Example as expl


V, choice_a, KL, w = expl.sample()

# Running a custom run with age-specific gridz ------------------------------

# Example gridspace...
grida_space = np.zeros(120)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 100
    elif a == 1:
        grida_space[a] = 200
    elif a < 5:
        grida_space[a] = grida_space[a-1]*1.5
    elif a < 10:
        grida_space[a] = grida_space[a-1]*1.4
    elif a < 20:
        grida_space[a] = grida_space[a-1]*1.2
    elif a < 40:
        grida_space[a] = grida_space[a-1]*1.1
    else:
        grida_space[a] = grida_space[a-1]*1.05

grida = np.concatenate((-np.flip(grida_space[0:(np.int(len(grida_space)/2))]),[0],grida_space))


gridz = values
Pi = transition
n = 40
xi = 1
delta = 0.05
beta = 0.98
alpha = 0.4
sigma_x = 1.74
sigma_y = 1.33
transfer = 13 * 89 #Simulating BF/unemployment/etc transfers

r = 0.05
w=1

# Running example with temptation
temptation = True
V, choice_a, KL, A = expl.age_specific(Pi, grida, gridz, temptation, \
                                       n, w, xi, r, delta, beta, alpha, sigma_x, sigma_y, transfer)

# Running example without temptation
temptation = False
nt_V, nt_choice_a, nt_KL, nt_A = expl.age_specific(Pi, grida, gridz, temptation, \
                                       n, w, xi, r, delta, beta, alpha, sigma_x, sigma_y, transfer)
#%% Menu

# Options here will include:
    
    # Setting the income process according to some inputable parameters:
        # Number of periods
        # Empirical data from where to base the process
        
    # Calibrating the model:
        # Externally calibrated preference parameters -> check resulting distribution
        # Internally calibrated preference parameters -> match wealth distribution / savings pattern
        
    # Simulating counterfactuals from established model
    
    # Plotting results
    
    # Exit program

#%% Calibration

import CalibrationB as cb

## ONLY POF ANALYSIS

# Loading empirical data

pof_path = "C:\Backup Acer\EESP-Mestrado\Dissertação\POF\R Project\Exports\\temptation_filtered.csv"
pof_df = cb.import_POF_data(pof_path)
data_x, data_y = cb.select_data(pof_df)

# Selecting a moment grid for calibration. Options: "arbitrary", "percentile5", "percentile10"
# gridm = cb.create_cons_moment_grid("percentile5")

# Selecting percentile grid for moment calibration
gridq = np.arange(0,1,0.1)

# Building which parameters will be explored:
grid_x0 = np.arange(60*12,90*12,5*12)
grid_sigx = np.arange(1.1,1.5,0.01)
grid_sigy = np.arange(1.01,1.3,0.01)
grid_xi = np.arange(1,5,0.1)

# Running the simulations
sqe = cb.run_x0_simulations(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, xi)
sol = np.argwhere(sqe == np.min(sqe))

# Running Two-Step simulations
sqe_I, sqe_W, W_new, Omega = cb.run_two_steps(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi)

#%% Simulating

import OLG_A_Power as olg

# Example gridspace...
grida_space = np.zeros(120)
for a in range(len(grida_space)):
    if a == 0: 
        grida_space[a] = 100
    elif a == 1:
        grida_space[a] = 200
    elif a < 5:
        grida_space[a] = grida_space[a-1]*1.5
    elif a < 10:
        grida_space[a] = grida_space[a-1]*1.4
    elif a < 20:
        grida_space[a] = grida_space[a-1]*1.2
    elif a < 40:
        grida_space[a] = grida_space[a-1]*1.1
    else:
        grida_space[a] = grida_space[a-1]*1.05

grida = np.concatenate((-np.flip(grida_space[0:(np.int(len(grida_space)/2))]),[0],grida_space))

gridz = values
Pi = transition
n = 40
xi = 0.452
delta = 0.05
beta = 0.98
alpha = 0.4
sigma_x = 1.6
sigma_y = 1.42
transfer = 0
x0 = 75*12

A = 1 # THINK ABOUT HOW TO DEAL WITH THIS
r_low = 0
r_high = 0.2
# Check if using correct mass grid for z
mass_z = np.concatenate((np.repeat(0.02,5),np.repeat(0.1,8),np.repeat(0.02,5))) 

# Computing a single GE - with Temptation
KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, r = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0, temptation = True, tol = 1e-2, maxit = 10000)

# Computing a single GE - without Temptation
KL_nt, w_nt, C_nt, x_nt, y_nt, V_nt, choice_a_nt, distr_mass_nt, k_mass_nt, k_total_nt, c_mass_nt, c_total_nt, r_nt = olg.general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0, temptation = False, tol = 1e-2, maxit = 10000)


#%% Plotting

import Graphics as grp

mass_by_k, mass_by_age_k = grp.capital_distributions(n, grida, gridz, distr_mass_nt, k_mass_nt)
grp.plot_total_k_distr(mass_by_k, grida, 10000, "Calibration A, With Temptation")
grp.plot_total_k_lorenz(mass_by_k, grida, "Calibration A, With Temptation")

# Plotting comparison: Temptation vs Non-Temptation
mass_by_k, mass_by_age_k = grp.capital_distributions(n, grida, gridz, distr_mass, k_mass)
mass_by_k_nt, mass_by_age_k_nt = grp.capital_distributions(n, grida, gridz, distr_mass_nt, k_mass_nt)

grp.compare_total_k_distr(mass_by_k, mass_by_k_nt, grida, bin_size = 2500, description = "Calibration (B)", label1 = "With Temptation", label2 = "Without Temptation", log = False, trim_upper = True, trim_value = 600000)
grp.compare_total_k_lorenz(mass_by_k, mass_by_k_nt, grida, description = "Calibration (B)", label1  = "With Temptation", label2 = "Without Temptation")

# Plotting Asset evolution by quantile
age_start = 25
quants = [0.2,0.5,0.9,0.95,0.99]

grp.plot_k_evolution(age_start, n, mass_by_age_k, quants, grida, description = "with temptation")
grp.plot_k_evolution(age_start, n, mass_by_age_k_nt, quants, grida, description = "without temptation")

#%% Calculating some stats

olg.calculate_wealth_gini(n, grida, gridz, distr_mass, k_mass)
olg.calculate_wealth_gini(n, grida, gridz, distr_mass_nt, k_mass_nt)

