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
#                            Temptation                                      #
##############################################################################

# Solves for allocation of between regular and temptation consumption

import numpy as np
from numba import njit

@njit
def calculate_U(x,sigma_x, x0, x_bound):
    
    if x > x_bound:
        U = (x-x0)**(1 - sigma_x) / (1-sigma_x)
    elif x > 0:
        U = (x_bound-x0)**(1 - sigma_x) / (1-sigma_x)
    else:
        U = -np.Inf
    return U

@njit
def calculate_T(y, sigma_y, xi, y_bound):
    
    if y > y_bound:
        T = xi*y**(1-sigma_y) / (1-sigma_y)
    elif xi == 0:
        T = 0
    else:
        T = xi*y_bound**(1-sigma_y) / (1-sigma_y)
    return T

@njit
def f(x, c, xi, sigma_x, sigma_y, x0):
    return x + (x-x0)**(sigma_x/sigma_y)*xi**(1/sigma_y) - c

@njit
def df(x, c, xi, sigma_x, sigma_y, x0):
    return 1 + xi**(1/sigma_y)*(sigma_x/sigma_y)*(x-x0)**((sigma_x-sigma_y)/sigma_y)

@njit
def Newton(f, df, x_init, e, maxit, c, xi, sigma_x, sigma_y, x0):
    xn = x_init
    for n in range(0,maxit):
        fxn = f(xn, c, xi, sigma_x, sigma_y, x0)
        if abs(fxn)<e:
            return xn
        dfxn = df(xn, c, xi, sigma_x, sigma_y, x0)
        if dfxn == 0:
            print("Reached zero derivative in iteration",n,"please review initial guess")
            return None
        xn = xn - fxn/dfxn
    print("Maximum iterations (",n,") reached. \n Current x:",xn,"\n Current f(x):",fxn)
    print("\n C evaluated: ", c,"\n sx evaluated: ", sigma_x,"\n sy evaluated: ", sigma_y)
    return None  

@njit
def create_allocation_grid_by_age(grida, gridz, xi, sigma_x, sigma_y, w, r, n, transfer, x0):
    
    print("Starting allocation grids by age...")
    
    # Starting timer
    # start_time = time.time()
    
    # Defining consumption spaces C, x, y (productivity, current asset, next asset, period)
    C = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    x = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    y = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    
    for a1 in range(len(grida)):
        for age in range(n):
            for z in range(len(gridz)):
                for a2 in range(len(grida)):
                    C[z,a1,a2,age] = w*gridz[z,age] + (1+r)*grida[a1] - grida[a2] + (transfer - w*gridz[z,age])*(w*gridz[z,age] < transfer)
                    if C[z,a1,a2,age] <= 1 or xi == 0:
                        C[z,a1,a2:,age] = max([C[z,a1,a2,age],0])
                        x[z,a1,a2:,age] = max([C[z,a1,a2,age],0])
                        y[z,a1,a2:,age] = 0
                        next
                    elif C[z,a1,a2,age] <= x0 + 1:
                        C[z,a1,a2:,age] = C[z,a1,a2,age]
                        x[z,a1,a2:,age] = C[z,a1,a2,age]
                        y[z,a1,a2:,age] = 0
                    else:
                        x[z,a1,a2,age] = Newton(f,df,(x0 + 0.0001),1e-6,int(1e6),C[z,a1,a2,age], xi, sigma_x, sigma_y, x0)
                        y[z,a1,a2,age] = C[z,a1,a2,age] - x[z,a1,a2,age]

        
        # Reporting evolution of calculation
        # print(np.round(a1/len(grida)*100,2),"% concluded \n")
    
    # hour = round(((time.time() - start_time)/60)//60)
    # minute = round((time.time() - start_time)//60 - hour*60)
    # second = round((time.time() - start_time) - hour*60*60 - minute*60)
                   
   #  print("Allocations calculated! \n Time elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")
    
    return C, x, y

@njit
def calculate_allocation_by_C(C, x0, sigma_x, sigma_y, xi):
    
    x = np.zeros(shape = len(C))
    y = np.zeros(shape = len(C))
    frac = np.zeros(shape = len(C))
    dxdc = np.zeros(shape = len(C))
    dxdc[0] = 1
    
    for i in range(len(C)):
        if C[i] <= 1 or xi == 0:
            C[i] = max([C[i],0])
            x[i] = max([C[i],0])
            y[i] = 0
            if i > 0:
                dxdc[i] = (x[i]-x[i-1]) / (C[i]-C[i-1])
            next
        elif C[i] <= x0 + 1:
            C[i] = C[i]
            x[i] = C[i]
            y[i] = 0
            if i > 0:
                dxdc[i] = (x[i]-x[i-1]) / (C[i]-C[i-1])
        else:
            x[i] = Newton(f,df,(x0 + 0.0001),1e-6,int(1e6),C[i], xi, sigma_x, sigma_y, x0)
            y[i] = C[i] - x[i]
            frac[i] = y[i] / C[i]
            if i > 0:
                dxdc[i] = (x[i]-x[i-1]) / (C[i]-C[i-1])
    
    return x, y, frac, dxdc





