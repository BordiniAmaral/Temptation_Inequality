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
def calculate_U(x,sigma_x, x0, bound):
    
    if x > x0 + bound:
        U = (x-x0)**(1 - sigma_x) / (1-sigma_x)
    elif x > 0:
        U = (bound)**(1 - sigma_x) / (1-sigma_x)
    else:
        U = -np.Inf
    return U

@njit
def calculate_T(y, sigma_y, xi, bound):
    
    if y > bound:
        T = xi*y**(1-sigma_y) / (1-sigma_y)
    else:
        T = xi*bound**(1-sigma_y) / (1-sigma_y)
    return T

'''
@njit
def y_allocation(y, *args):
    c, xi, sigma_x, sigma_y = args
    residual = c - y - (xi * ((sigma_y-1)/sigma_y) * (sigma_x/(sigma_x-1)) * y**(-1/sigma_y))**(-sigma_x)
    return residual
'''

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

'''
def create_allocation_grid(grida,gridz,xi,sigma_x,sigma_y, w, r):
    
    # Starting timer
    start_time = time.time()
    
    # Defining consumption spaces C,x,y(productivity, current asset, next asset)
    C = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    x = np.zeros(shape = (len(gridz), len(grida), len(grida)))
    y = np.zeros(shape = (len(gridz), len(grida), len(grida)))

    # Defining last period consumption space C,x,y_last(current asset)
    C_last = np.zeros(shape = (len(grida)))
    x_last = np.zeros(shape = (len(grida)))
    y_last = np.zeros(shape = (len(grida)))
    
    for a1 in range(len(grida)):
        C_last[a1] = (1+r)*grida[a1]
        if C_last[a1]>0:
            y_last[a1] = fsolve(func = y_allocation, x0 = C_last[a1]/2, args = (C_last[a1],xi,sigma_x,sigma_y))
            x_last[a1] = C_last[a1] - y_last[a1]
        else:
            y_last[a1] = 0
            x_last[a1] = 0
            
        for z in range(len(gridz)):
            still_in_budget = True
            for a2 in range(len(grida)):
                if still_in_budget:
                    C[z,a1,a2] = w*gridz[z] + (1+r)*grida[a1] - grida[a2]
                    if C[z,a1,a2] < 0:
                        still_in_budget = False
                        C[z,a1,a2] = 0
                        x[z,a1,a2] = 0
                        y[z,a1,a2] = 0
                    else:
                        y[z,a1,a2] = fsolve(func = y_allocation, x0 = C[z,a1,a2]/2, args = (C[z,a1,a2],xi,sigma_x,sigma_y))
                        x[z,a1,a2] = C[z,a1,a2] - y[z,a1,a2]
                else:
                    C[z,a1,a2] = 0
                    x[z,a1,a2] = 0
                    y[z,a1,a2] = 0
    
    hour = round(((time.time() - start_time)/60)//60)
    minute = round((time.time() - start_time)//60 - hour*60)
    second = round((time.time() - start_time) - hour*60*60 - minute*60)
                   
    print(f"Allocations calculated! \n Time elapsed (Total): {hour}h {minute}min {second}s \n \n")
    
    return x, y, x_last, y_last
'''    

@njit
def create_allocation_grid_by_age(grida, gridz, xi, sigma_x, sigma_y, w, r, n, transfer, x0):
    
    print("Starting allocation grids by age...")
    
    # Starting timer
    # start_time = time.time()
    
    # Defining consumption spaces C,x,y(productivity, current asset, next asset)
    C = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    x = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    y = np.zeros(shape = (len(gridz), len(grida), len(grida), n))
    
    for a1 in range(len(grida)):
        for age in range(n):
            for z in range(len(gridz)):
                for a2 in range(len(grida)):
                    C[z,a1,a2,age] = w*gridz[z,age] + (1+r)*grida[a1] - grida[a2] + (transfer - w*gridz[z,age])*(w*gridz[z,age] < transfer)
                    if C[z,a1,a2,age] <= 1:
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
                        
                        # y[z,a1,a2,age] = fsolve(func = y_allocation, x0 = C[z,a1,a2,age]/10, args = (C[z,a1,a2,age],xi,sigma_x,sigma_y))
                        # x[z,a1,a2,age] = C[z,a1,a2,age] - y[z,a1,a2,age]
        
        # Reporting evolution of calculation
        # print(np.round(a1/len(grida)*100,2),"% concluded \n")
    
    # hour = round(((time.time() - start_time)/60)//60)
    # minute = round((time.time() - start_time)//60 - hour*60)
    # second = round((time.time() - start_time) - hour*60*60 - minute*60)
                   
   #  print("Allocations calculated! \n Time elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")
    
    return C, x, y


########################## TESTING NUMBA ##################################

# start = time.time()
# x, y = create_allocation_grid_by_age(grida,gridz,xi,sigma_x,sigma_y, w, r, n)
# end = time.time()
# print("Elapsed (with compilation) = %s" % (end - start))

# start = time.time()
# Newton(f,df,1,1e-10,int(1e6),3000, xi, sigma_x, sigma_y)
# end = time.time()
# print("Elapsed (with compilation) = %s" % (end - start))





