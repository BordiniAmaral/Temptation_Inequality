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
#                            OLG Solution                                    #
##############################################################################

# Takes parameters as given and returns elements of interest:
    # Distribution per cohort
    # Policy functions

import numpy as np
from numba import njit
import Temptation_CRRA as tpt


@njit
def calculate_KL(r, delta, alpha, A):
    return ((r+delta)/(alpha*A))**(1/(alpha-1))

@njit
def calculate_w(alpha, A, KL):
    return (1-alpha)*A*KL**alpha

# Computes aggregate Labor measure
@njit
def calculate_labor(gridz, mass_z, n):
    
    L_mass = np.zeros(shape = (n,len(gridz)))
    
    for z in range(len(gridz)):
        current_z_distr = mass_z[z]
        for age in range(n):
            L_mass[age,z] = current_z_distr * gridz[z,age] * (1/n)
    
    return np.sum(L_mass)

@njit
def compute_utility_grids(x,y,grida,gridz,n,sigma_x,sigma_y,xi,x0, bound):
    
    # Defining regular utility and temptation utility space U(x) and T(y)
    U = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    T = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    
    for z in range(len(gridz)):
        for a1 in range(len(grida)):
            for a2 in range(len(grida)):
                for age in range(n):
                    U[z,a1,a2,age] = tpt.calculate_U(x[z,a1,a2,age], sigma_x, x0, bound)
                    T[z,a1,a2,age] = tpt.calculate_T(y[z,a1,a2,age], sigma_y, xi, bound)
    
    return U, T

@njit
def partial_sol_age_specific(n, beta, alpha, Pi, gridz, grida, x, y, sigma_x, sigma_y, xi, x0, temptation):
    
    print("Computing partial equilibrium age-specific solution...")
    # Defining value function space V(period, productivity, current asset)
    V = np.zeros(shape = (n, len(gridz), len(grida)))
    
    bound = 1 # Bounding CRRA by evaluating only up to a certain distance from minimum value
    U_bound = (bound)**(1 - sigma_x) / (1-sigma_x)
    
    print("Computing utility grids...")
    U, T = compute_utility_grids(x, y, grida ,gridz, n, sigma_x, sigma_y, xi, x0, bound)

    # Creating choice objects choice (period, productivity, current asset)
    choice_a = np.zeros(shape = (n, len(gridz), len(grida)))
    
    # Creating options (productivity, current asset, next asset, age)
    # from which to draw the maximum
    Aux = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    
    # Creating feasible a2 indicator
    # Allows for adding borrowing constraints, if desired!
    possible = U > -np.Inf
    
    zero_asset_index = np.int64(np.argwhere(grida == 0)[0][0])
    
    print("Starting backards iteration of optimal solutions")
    for period in range(n-1, -1, -1):
        
        # Last period implies consumption of all positive assets
        if period == n-1:
            for a1 in range(len(grida)):
                choice_a[period, :, a1] = zero_asset_index
                for z1 in range(len(gridz)):
                    if possible[z1,a1,zero_asset_index,period]:
                        V[period, z1, a1] = U[z1,a1,zero_asset_index, period] + T[z1,a1,zero_asset_index, period]
                    else:
                        V[period, z1, a1] = -np.Inf
                        
        # Other periods occur normally
        else:
            
            # Computing V_next(z1,a2) separately for an improved loop in what
            # follows
            V_next = np.zeros(shape = (len(gridz),len(grida)))
            for z1 in range(len(gridz)):
                for a2 in range(len(grida)):
                    if V[period+1, 0, a2] == -np.Inf:
                        V_next[z1,a2] = -np.Inf
                    else:
                        for z2 in range(len(gridz)):
                            V_next[z1,a2] = V_next[z1,a2] + Pi[z1,z2]*(V[period+1,z2,a2] - T[z2,a2,int(choice_a[period+1,z2,a2]),period]*temptation)
            
            Aux = np.zeros(shape = (len(gridz), len(grida), len(grida)))
            # Trying to cover the state space in a smart way, avoiding
            # computing regions that are clearly non optimal
            for a1 in range(len(grida)):
                for z1 in range(len(gridz)):
                    
                    # Excluding impossible steps straightforwardly
                    a_possible_index = np.sum(possible[z1,a1,:,period])
                    Aux[z1,a1,a_possible_index:] = -np.Inf
                    
                    minimum_indexes = np.argwhere(U[z1,a1,:,period] == U_bound)
                    unfeasible_indexes = np.argwhere(U[z1,a1,:,period] == -np.Inf)
                    
                    # Starting loop on a2 from highest possible, then detecting
                    # when Aux starts to move away from (local) maximum.
                    # Rule (arbitrary): 5 consecutive decreases to Aux with
                    # decreasing a2 below minimum consumption
                    decreasing = 0
                    for a2 in range(a_possible_index-1,-1,-1): ##  range(len(grida))
                        Aux[z1,a1,a2] = U[z1,a1,a2,period] + T[z1,a1,a2,period] + beta*V_next[z1,a2]
                        
                        if a2 < a_possible_index-1:
                            if (minimum_indexes.size == 0):
                                if Aux[z1,a1,a2] < Aux[z1,a1,a2+1]:
                                    decreasing = decreasing + 1
                                else:
                                    decreasing = 0
                            elif a2 < minimum_indexes[0][0]:
                                if Aux[z1,a1,a2] < Aux[z1,a1,a2+1]:
                                    decreasing = decreasing + 1
                                else:
                                    decreasing = 0
                        
                        if decreasing >= 5:
                            Aux[z1,a1,:a2] = -np.Inf
                            break
                    
                    # Applying rule:
                        # For positive future assets, picking only the highest feasible option
                        # For negative future assets, picking only the the feasible options below zero_asset (i.e. the agent pays up before being able to consume)
                    
                    #print("FOR DEBUG: (z1,a1,period) = ",z1,a1,period)
                    
                    if minimum_indexes.size != 0 and unfeasible_indexes.size != 0:
                        cut_index = max([(minimum_indexes[0][0]+1), min([zero_asset_index + 1,(unfeasible_indexes[0][0])])])
                        Aux[z1,a1,cut_index:] = -np.Inf
                    elif minimum_indexes.size != 0 and unfeasible_indexes.size == 0:
                        cut_index = minimum_indexes[0][0]+1
                        Aux[z1,a1,cut_index:] = -np.Inf
                    elif minimum_indexes.size == 0 and unfeasible_indexes.size != 0: 
                        cut_index = unfeasible_indexes[0][0]
                        Aux[z1,a1,cut_index:] = -np.Inf
                    
                    V[period,z1,a1] = max(Aux[z1,a1,:])
                    choice_a[period,z1,a1] = np.argmax(Aux[z1,a1,:])

    print("Backwards iteration concluded ")
    return V, choice_a

@njit
def aggregate_capital(grida, distr_mass, gridz, n):
    
    k_mass = np.zeros(shape = (n+1, len(gridz), len(grida)))
    
    for a in range(len(grida)):
        current_a = grida[a]
        for age in range(n+1):
            for z in range(len(gridz)):
                k_mass[age,z,a] = distr_mass[age,z,a] * current_a
    
    return k_mass, np.sum(k_mass)

@njit
def aggregate_consumption(distr_mass, grida, gridz, n, C, choice_a):
    
    c_mass = np.zeros(shape = (n, len(gridz), len(grida)))
    
    for a in range(len(grida)):
        for age in range(n):
            for z in range(len(gridz)):
                c_mass[age,z,a] = distr_mass[age,z,a] * C[z,a,np.int32(choice_a[age,z,a]),age]
    
    return c_mass, np.sum(c_mass)

@njit
def stationary_distribution(choice_a, grida, gridz, mass_z, Pi, n):
    
    print("Computing stationary distribution...")
    distr_mass = np.zeros(shape = (n+1, len(gridz), len(grida)))
    zero_asset_index = np.int32(np.argwhere(grida == 0)[0][0])
    
    for age in range(n+1):
        if age == 0:
            distr_mass[age,:,zero_asset_index] = mass_z
        else:
            for z1 in range(len(gridz)):
                for a1 in range(len(grida)):
                    current_a = np.int64(choice_a[age-1,z1,a1])
                    for z2 in range(len(gridz)):
                        distr_mass[age,z2,current_a] = distr_mass[age,z2,current_a] + Pi[z1,z2] * distr_mass[age-1,z1,a1]
    
    return distr_mass / (n+1)

# Takes as input parameters, runs simulations looking for GE interest rate
# Returns resulting distributions and equilibrium parameters
@njit
def general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, A, x0,
                        temptation, tol = 1e-2, maxit = 100):
    
    print("\nRunning initial calculations...")
    # Starting counter and error
    i = 0
    
    zsum = calculate_labor(gridz, mass_z, n)
    
    excess_KL = np.zeros(shape = (2))
    init_r = np.array([r_low, r_high])
    window = (init_r[1]-init_r[0])
    
    print("Checking given guesses...")
    # Running first for inital guesses
    
    adequate_r = False
    both_tested = False
    reduced = False
    
    # Testing lower r (and updating upper r when appropriate)
    while not adequate_r:
        print("\nTesting r:",init_r[0])
        KLd, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total = run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, init_r[0], mass_z, transfer, A, temptation, x0)
        
        L_total = zsum / w # zsum is total salary mass considering the ranges in POF
        KLs = k_total / L_total
        
        if k_total > 0:
            excess_KL[0] = KLs - KLd
            if excess_KL[0] < 0:
                print("\nLower r checked. \nValue:",init_r[0],"\nCurrent error:",excess_KL[0])
                adequate_r = True
                if excess_KL[1] > 0:
                    both_tested = True
            else:
                excess_KL[0] = 0
                excess_KL[1] = KLs - KLd
                init_r = np.array([init_r[0] - window/2, init_r[0]])
                print("\nLower r too high.\nReducing range by half-window to: ", np.round(init_r[0],2),"-", np.round(init_r[1],2))
                reduced = True
        else:
            if reduced:
                print("\nReduced lower r guess too much. Reducing window before next increase.")
                window = window/2
            init_r = np.array([init_r[0] + window/2, init_r[1] + window/2])
            print("\nNegative capital supply reached.\nIncreasing range by half-window to:", np.round(init_r[0],2),"-", np.round(init_r[1],2))
    
    # Testing upper r, if not done in first part
    if not both_tested:
        adequate_r = False
        while not adequate_r:
            print("\nTesting r:",init_r[1])
            KLd, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total = run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, init_r[1], mass_z, transfer, A, temptation, x0)
            
            L_total = zsum / w # zsum is total salary mass considering the ranges in POF
            KLs = k_total / L_total
            
            if k_total > 0:
                excess_KL[1] = KLs - KLd
                if excess_KL[1] > 0:
                    print("\nUpper r checked. \nValue:",init_r[1],"\nCurrent error:",excess_KL[1])
                    adequate_r = True
                else:
                    excess_KL[1] = 0
                    excess_KL[0] = KLs - KLd
                    init_r = np.array([init_r[1], init_r[1]+ window/2])
                    print("\nUpper r too low.\nIncreasing range by half-window to: ", np.round(init_r[0],2),"-", np.round(init_r[1],2))
    
    print("\n----------------------------\nInitial guesses finished.\nRange:", np.round(init_r[0],2),"-", np.round(init_r[1],2))
    
    r_low = init_r[0]
    r_high = init_r[1]
    
    if (excess_KL[0] < 0) and (excess_KL[1] > 0):
        print("\n--------------------------------------\n        Proceeding to Bissection \n--------------------------------------")
        excess_KL = 1
        
        while (np.abs(excess_KL) > tol) and (i < maxit):
            
            print("Starting iteration",i+1,"\n")
            
            r = (r_high + r_low)/2 
         
            KLd, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total = run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r, mass_z, transfer, A, temptation, x0)
            
            L_total = zsum / w
            KLs = k_total / L_total
            
            excess_KL = KLs - KLd
            i = i+1
            
            print("\nIteration",i,"concluded.\n Current r interval:",r_low,"and",r_high,"\n Current r guess:",r,"\n Current error:",excess_KL,"\n--------------------------------------------\n")
            
            # Updating bissection interval
            if excess_KL < 0:
                r_low = r
            else:
                r_high = r
            
        
        if i == maxit:
            print("Maximum iterations reached")
        
        if np.abs(excess_KL) < tol:
            print("Equilibrium calculated. Values:\n Equilibrium r:", r)
    else:
        return None
    
    return KLd, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, r, init_r

@njit
def run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r, mass_z, transfer, A, temptation, x0):
    
    KL = calculate_KL(r, delta, alpha, A)
    w = calculate_w(alpha, A, KL)
    
    # Adapting gridz to correspond to hh annual disposable income in R$:
    gridz_new = gridz/w
    
    # Creating allocation grid
    C, x, y = tpt.create_allocation_grid_by_age(grida, gridz_new, xi, sigma_x, sigma_y, w, r, n, transfer, x0)
    
    # Solving lifecycle
    V, choice_a = partial_sol_age_specific(n, beta, alpha, Pi, gridz, grida, x, y, sigma_x, sigma_y, xi, x0, temptation)
    
    # Aggregating capital and consumption
    distr_mass = stationary_distribution(choice_a, grida, gridz_new, mass_z, Pi, n)
    k_mass, k_total = aggregate_capital(grida, distr_mass, gridz_new, n)
    c_mass, c_total = aggregate_consumption(distr_mass, grida, gridz_new, n, C, choice_a)
    
    return KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total

def calculate_wealth_gini(n, grida, gridz, distr_mass, k_mass):
    
    mass_by_k = np.zeros(shape = len(grida))
    mass_by_age_k = np.zeros(shape = (n, len(grida)))
    
    for a in range(len(grida)):
        mass_by_k[a] = np.sum(distr_mass[:,:,a])
        for n in range(n):
            mass_by_age_k[n,a] = np.sum(distr_mass[n,:,a])
            
    cum_mass_by_k = np.cumsum(mass_by_k)
    cum_kfrac_by_k = np.cumsum(mass_by_k * grida) / np.sum(mass_by_k * grida)
    
    cum_mass_by_k_lag = np.concatenate(([0],cum_mass_by_k[:-1]))
    
    gini = 1 - np.sum((cum_mass_by_k - cum_mass_by_k_lag)*cum_kfrac_by_k) / 0.5
    
    print("Gini index calculated:",np.round(gini,decimals=3))
    
    return gini

# Running General Equilibrium results with different betas in order to match
# a certain level of aggregate capital

@njit
def ge_match_by_capital(beta0, step, tol, k_target, n, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r, mass_z, transfer, A, x0, temptation):
    
    # Here, we consider beta0 the one used without temptation, and plain CRRA. Thus, it is an
    # lower bound to the beta we are looking for. First I shall test step-incremented
    # {1,...,m} higher betas until I find one beta_{m} s.t. equilibrium aggregate capital
    # with given known r is lower than the target. Then, I shall execute a bissection 
    # with the range beta_{m-1} - beta{m} until the relative error is below tolerance
    
    # Finding the beta_m and conducing bissection
    m = 1
    error = 1
    bissection = False
    count = 1
    beta_low = beta0
    beta_high = beta0 + m*step
    
    while abs(error) > tol:
        
        if not bissection:
            beta = beta0 + m*step
        else:
            beta = (beta_high + beta_low)/2
        print("\n** Starting iteration ",count,"of capital matching**\n     Current guess:",beta)
        
        KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total = run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r, mass_z, transfer, A, temptation, x0)
        
        error = (k_total - k_target) / k_target
        
        if error > 0 and not bissection:
            print("\nUpper bound for beta found: ",np.round(beta,3),". Error = ",np.round(error,6),".\nProceeding to bissection.")
            bissection = True
            beta_high = beta0+step*m
            beta_low = beta0+step*(m-1)
        elif error < 0 and not bissection:
            m = m+1
            print("\nGE capital still below target. Error = ",np.round(error,6),".\nIncreasing beta to: ",np.round(beta0+step*m,3))
        elif abs(error) > tol and bissection:
            print("\Iteration concluded.\nCurrent beta: ",np.round(beta,3),"\nCurrent error: ",np.round(error,6))
            if error > 0:
                print("\nDecreasing beta.")
                beta_high = beta
            else:
                print("\nIncreasing beta.")
                beta_low = beta
        elif abs(error) < tol and bissection:
            print("\n\n**Bissection concluded**\nBeta = ",np.round(beta,3),"\nError = ",np.round(error,6))
        count = count + 1
    
    results = (KL, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, r)
    
    return beta, results
            
def compute_k_supply_curve(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, grid_r, r_ref, mass_z, transfer, A, temptation, x0):
    
    k_supply = np.zeros(len(grid_r))
    
    # Computing reference GE labor supply
    KL = calculate_KL(r_ref, delta, alpha, A)
    w = calculate_w(alpha, A, KL)
    
    # Adapting gridz to correspond to hh annual disposable income in R$ in reference scenario:
    gridz_new = gridz/w
    
    for r in range(len(grid_r)):
        
        # Creating allocation grid
        C, x, y = tpt.create_allocation_grid_by_age(grida, gridz_new, xi, sigma_x, sigma_y, w, grid_r[r], n, transfer, x0)
        
        # Solving lifecycle
        V, choice_a = partial_sol_age_specific(n, beta, alpha, Pi, gridz, grida, x, y, sigma_x, sigma_y, xi, x0, temptation)
        
        # Aggregating capital and consumption
        distr_mass = stationary_distribution(choice_a, grida, gridz_new, mass_z, Pi, n)
        k_mass, k_total = aggregate_capital(grida, distr_mass, gridz_new, n)
        
        k_supply[r] = k_total
        
    return k_supply

# Remember below that A = 1/w by construction of our calibration
def compute_k_demand_curve(grid_r, delta, w, A, gridz, mass_z, n, alpha):
    
    k_demand = np.zeros(len(grid_r))
    
    zsum = calculate_labor(gridz, mass_z, n)
    L_total = zsum / w
    
    for r in range(len(grid_r)):
        
        k_demand[r] = ((grid_r[r]+delta)/(alpha*A))**(1/(alpha-1))*L_total
    
    return k_demand