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
#                            Sample Run                                      #
##############################################################################


# Provides sample run for testing the functions

def sample():
    
    import Temptation as tpt
    import OLG as olg
    import numpy as np
    
    # Example parameters
    n = 40
    xi = 1
    r = 0.05
    delta = 0.05
    beta = 0.98
    alpha = 0.4
    sigma_x = 1.74
    sigma_y = 1.33
    
    # Example gridspace
    grida = np.linspace(-20,50,71)
    gridz = np.linspace(0.5,5,10)
    
    # Example income Markov chain: P(z_column | z_line)
    Pi = np.array([[0.5,0.35,0.15,0,0,0,0,0,0,0],
                   [0.2,0.4,0.3,0.1,0,0,0,0,0,0],
                   [0.05,0.2,0.5,0.2,0.05,0,0,0,0,0],
                   [0,0.05,0.2,0.5,0.2,0.05,0,0,0,0],
                   [0,0,0.05,0.2,0.5,0.2,0.05,0,0,0],
                   [0,0,0,0.05,0.2,0.5,0.2,0.05,0,0],
                   [0,0,0,0,0.05,0.2,0.5,0.2,0.05,0],
                   [0,0,0,0,0,0.05,0.2,0.5,0.2,0.05],
                   [0,0,0,0,0,0,0.1,0.3,0.4,0.2],
                   [0,0,0,0,0,0,0,0.15,0.35,0.5]])
    
    KL = olg.compute_KL(r, delta, alpha)
    w = olg.compute_wage(KL, alpha)
    x, y, x_last, y_last = tpt.create_allocation_grid(grida,gridz,xi,sigma_x,sigma_y, w, r)
    V, choice_a = olg.partial_sol_baseline(n, beta, alpha, Pi, gridz, grida, x, y, x_last, y_last, sigma_x, sigma_y, xi, temptation = True)
    
    return V, choice_a, KL, w

def age_specific(Pi,
                grida,
                gridz,
                temptation,
                n = 40,
                w=1,
                xi = 1,
                r = 0.05,
                delta = 0.05,
                beta = 0.98,
                alpha = 0.4,
                sigma_x = 1.74,
                sigma_y = 1.33,
                transfer = 0):
    
    import Temptation as tpt
    import OLG as olg
    
    A = olg.compute_A(r, delta, alpha)
    KL = olg.compute_KL(r, delta, alpha, A)
    x, y = tpt.create_allocation_grid_by_age(grida,gridz,xi,sigma_x,sigma_y, w, r, n, transfer)
    V, choice_a = olg.partial_sol_age_specific(n, beta, alpha, Pi, gridz, grida, x, y, sigma_x, sigma_y, xi, temptation)
    
    return V, choice_a, KL, A
