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
import Temptation as tpt

# @njit
# def calculate_KL(r, delta, alpha, A):
#     return ((r+delta)/(alpha*A))**(1/alpha-1)

@njit
def calculate_A(r, delta, alpha, KL):
    return ((r+delta)/alpha)*(KL)**(1-alpha)

@njit
def calculate_w(alpha, A, KL):
    return (1-alpha)*A*KL**alpha

@njit
def calculate_labor(gridz, mass_z, n):
    
    L_mass = np.zeros(shape = (n,len(gridz)))
    
    for z in range(len(gridz)):
        current_z_distr = mass_z[z]
        for age in range(n):
            L_mass[age,z] = current_z_distr * gridz[z,age] * (1/n)
    
    return np.sum(L_mass)

@njit
def partial_sol_baseline(n, beta, alpha, Pi, gridz, grida, x, y, x_last, y_last, sigma_x, sigma_y, xi, x0, temptation = True):
    
    print("\n Computing partial equilibrium baseline solution \n \n")
    # Defining value function space V(period, productivity, current asset)
    V = np.zeros(shape = (n, len(gridz), len(grida)))
    
    # Defining regular utility and temptation utility space U(x) and T(y)
    U = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    T = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    
    for z in range(len(gridz)):
        for a1 in range(len(grida)):
            for a2 in range(len(grida)):
                U[z,a1,a2] = calculate_U(x[z,a1,a2], sigma_x, x0)
                T[z,a1,a2] = calculate_T(y[z,a1,a2], sigma_y, xi)
    
    # Defining last period utilities 
    U_last = np.zeros(shape = (len(grida)))
    T_last = np.zeros(shape = (len(grida)))
    
    for a in range(len(grida)):
        U_last[a] = calculate_U(x_last[a], sigma_x, x0)
        T_last[a] = calculate_T(y_last[a], sigma_y, xi)
    
    # Creating choice objects choice(period, productivity, current asset)
    choice_a = np.zeros(shape = (n, len(gridz), len(grida)))
    
    # Creating options (productivity, current asset, next asset)
    # from which to draw the maximum
    Aux = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    
    # Creating feasible a2 indicator
    # Allows for adding borrowing constraints, if desired!
    possible = U > 0
    possible_last = U_last > 0
    
    for period in range(n-1, -1, -1):
        
        # Last period has no labor
        if period == n-1:
            for a in range(len(grida)):
                V[period, :, a] = U_last[a] + T_last[a]
                choice_a[period, :, a] = np.argwhere(grida == 0)
        
        # Previous-than-last period sees the next as the final
        elif period == n-2:
            for a2 in range(len(grida)):
                if not(possible_last[a2]):
                    Aux[:,:,a2] = 0
                else:
                    V_next = V[period+1,1,a2] - T_last[a2]*temptation
                    for z1 in range(len(gridz)):
                        for a1 in range(len(grida)):
                            if possible[z1,a1,a2]:
                                Aux[z1,a1,a2] = U[z1,a1,a2] + T[z1,a1,a2] + beta*V_next
                            else:
                                Aux[z1,a1,a2] = 0
            
            for z1 in range(len(gridz)):
                for a1 in range(len(grida)):
                    V[period,z1,a1] = max(Aux[z1,a1,:])
                    choice_a[period,z1,a1] = np.argmax(Aux[z1,a1,:])
                        
        # Other periods occur normally
        else:
            for a2 in range(len(grida)):
                for z1 in range(len(gridz)):
                    if min(V[period+1,Pi[z1,:]>0,a2]) == 0: #checking next period's feasibility
                        Aux[z1,:,a2] = 0
                    else:
                        V_next = 0
                        for z2 in range(len(gridz)):
                            V_next = V_next + Pi[z1,z2]*(V[period+1,z2,a2] - T[z2,a2,int(choice_a[period+1,z2,a2])]*temptation)
                        for a1 in range(len(grida)):
                            if possible[z1,a1,a2]:
                                Aux[z1,a1,a2] = U[z1,a1,a2] + T[z1,a1,a2] + beta*V_next
                            else:
                                Aux[z1,a1,a2] = 0
            
            for z1 in range(len(gridz)):
                for a1 in range(len(grida)):
                    V[period,z1,a1] = max(Aux[z1,a1,:])
                    choice_a[period,z1,a1] = np.argmax(Aux[z1,a1,:])
    
    return V, choice_a

@njit
def compute_utility_grids(x,y,grida,gridz,n,sigma_x,sigma_y,xi,x0):
    
    # Defining regular utility and temptation utility space U(x) and T(y)
    U = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    T = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    
    for z in range(len(gridz)):
        for a1 in range(len(grida)):
            for a2 in range(len(grida)):
                for age in range(n):
                    U[z,a1,a2,age] = calculate_U(x[z,a1,a2,age], sigma_x, x0)
                    T[z,a1,a2,age] = calculate_T(y[z,a1,a2,age], sigma_y, xi)
    
    return U, T

@njit
def partial_sol_age_specific(n, beta, alpha, Pi, gridz, grida, x, y, sigma_x, sigma_y, xi, x0, temptation = True):
    
    print("Computing partial equilibrium age-specific solution...")
    # Defining value function space V(period, productivity, current asset)
    V = np.zeros(shape = (n, len(gridz), len(grida)))
    
    print("Computing utility grids...")
    U, T = compute_utility_grids(x,y,grida,gridz,n,sigma_x,sigma_y,xi , x0)

    # Creating choice objects choice(period, productivity, current asset)
    choice_a = np.zeros(shape = (n, len(gridz), len(grida)))
    
    # Creating options (productivity, current asset, next asset, age)
    # from which to draw the maximum
    Aux = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    
    # Creating feasible a2 indicator
    # Allows for adding borrowing constraints, if desired!
    possible = U > 0
    
    zero_asset_index = np.int64(np.argwhere(grida == 0)[0][0])
    
    print("Starting backards iteration of optimal solutions")
    for period in range(n-1, -1, -1):
        
        # Last period implies consumption of positive assets
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
            for a2 in range(len(grida)):
                for z1 in range(len(gridz)):
                    V_next = 0
                    for z2 in range(len(gridz)):
                        if V[period+1, z2, a2] == -np.Inf:
                            Aux[z1,:,a2] = -np.Inf
                            break
                        else:
                            V_next = V_next + Pi[z1,z2]*(V[period+1,z2,a2] - T[z2,a2,int(choice_a[period+1,z2,a2]),period]*temptation)
                    
                    if V[period+1, 0, a2] == -np.Inf:
                        break
                    else:
                        for a1 in range(len(grida)):
                            if not possible[z1,a1,a2,period]:
                                Aux[z1,a1,a2] = -np.Inf
                            else:
                                Aux[z1,a1,a2] = U[z1,a1,a2,period] + T[z1,a1,a2,period] + beta*V_next
            
            for z1 in range(len(gridz)):
                for a1 in range(len(grida)):
                    # Applying rule:
                        # For positive assets, picking only the highest feasible option
                        # For negative assets, picking only the the feasible options below zero_asset (i.e. the agent pays up before being able to consume)
                    
                    #print("FOR DEBUG: (z1,a1,period) = ",z1,a1,period)
                    minimum_indexes = np.argwhere(U[z1,a1,:,period] == 0.01)
                    unfeasible_indexes = np.argwhere(U[z1,a1,:,period] == 0)
                    
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
def general_equilibrium(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r_low, r_high, mass_z, transfer, KL, x0,
                        temptation = True, tol = 1e-2, maxit = 100):
    
    print("\n Running initial calculations...")
    # Starting counter and error
    i = 0
    
    excess_KL = np.zeros(shape = (2))
    init_r = np.array([r_low, r_high])
    
    print("Checking given guesses...")
    # Running first for 
    for r in range(len(init_r)):
        print("\n Testing r:",init_r[r])
        A, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, gridz_new = run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, init_r[r], mass_z, transfer, KL, temptation, x0)
        
        L_total = calculate_labor(gridz_new, mass_z, n)
        KLs = k_total / L_total
        
        if k_total > 0:
            excess_KL[r] = KLs - KL
        else:
            print("Negative capital supply reached. Reasses given r.")
            excess_KL[r] = 1e10
    
    if excess_KL[0] < 0:
        print("Lower r checked: OK. Current error:",excess_KL[0])
    else:
        print("Lower r inadequate. Current error:",excess_KL[0],". Pick lower r (you may adopt this as your upper)")
    
    if excess_KL[1] > 0:
        print("Upper r checked: OK . Current error:",excess_KL[1])
    else:
        print("Upper r inadequate. Current error:",excess_KL[1],". Pick higher r (you may adopt this as your lower)")
    
    if (excess_KL[0] < 0) and (excess_KL[1] > 0):
        print("\n--------------------------------------\n        Proceeding to Bissection \n--------------------------------------")
        excess_KL = 1
        
        while (np.abs(excess_KL) > tol) and (i < maxit):
            
            print("Starting iteration",i,"\n")
            
            r = (r_high + r_low)/2 
         
            A, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, gridz_new = run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r, mass_z, transfer, KL, temptation, x0)
            
            L_total = calculate_labor(gridz_new, mass_z, n)
            KLs = k_total / L_total
            
            excess_KL = KLs - KL
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
    
    return A, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, r

@njit
def run_once(n, beta, delta, alpha, Pi, gridz, grida, sigma_x, sigma_y, xi, r, mass_z, transfer, KL, temptation, x0):
    
    A = calculate_A(r, delta, alpha, KL)
    w = calculate_w(alpha, A, KL)
    
    # Adapting gridz to correspond to salary in R$:
    gridz_new = gridz/w
    
    # Creating allocation grid
    C, x, y = tpt.create_allocation_grid_by_age(grida, gridz_new, xi, sigma_x, sigma_y, w, r, n, transfer, x0)
    
    # Solving lifecycle
    V, choice_a = partial_sol_age_specific(n, beta, alpha, Pi, gridz_new, grida, x, y, sigma_x, sigma_y, xi, x0, temptation)
    
    # Aggregating capital and consumption
    distr_mass = stationary_distribution(choice_a, grida, gridz_new, mass_z, Pi, n)
    k_mass, k_total = aggregate_capital(grida, distr_mass, gridz_new, n)
    c_mass, c_total = aggregate_consumption(distr_mass, grida, gridz_new, n, C, choice_a)
    
    return A, w, C, x, y, V, choice_a, distr_mass, k_mass, k_total, c_mass, c_total, gridz_new