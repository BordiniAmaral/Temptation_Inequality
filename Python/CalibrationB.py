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
#                          Calibration                                       #
##############################################################################

from numba import njit
import numpy as np
import pandas as pd

def import_POF_data(file_path):
    
    
    pof_df = pd.read_csv(file_path, sep = ";")
    
    return pof_df

@njit
def run_simulations_by_x0(grid_sigx, grid_sigy, x0, grid_x0, data_x, data_y, gridm, grid_xi, consumpt, cons_bin, y_data_sel_sum, W):
    
    print("------------------------ \n   Current x0:",x0,"  \n------------------------\n")
    sqe = np.zeros(shape = (len(grid_sigx), len(grid_sigy), len(grid_xi)))
    counter = 0
    sigy_start = np.int32(0)
    x0_index = np.argwhere(grid_x0 == x0)[0][0]
    
    for xi in range(len(grid_xi)):
        sigy_start = np.int32(0)
        for sigx in range(len(grid_sigx)):
            for sigy in range(sigy_start, len(grid_sigy)):
                
                if grid_sigy[sigy] > grid_sigx[sigx]:
                    sqe[sigx,sigy:,xi] = np.Inf #arbitrary high number
                    break
                
                y_sim = np.zeros(shape = (len(data_x)))
                x_sim = np.zeros(shape = (len(data_x)))
                
                for x in range(len(data_x)):
                    if consumpt[x] < (x0 + 2):
                        x_sim[x] = consumpt[x]
                    else:
                        x_sim[x] = Newton(f, df, (x0 + 0.0001), 1e-6,int(1e6), consumpt[x], grid_xi[xi], grid_sigx[sigx], grid_sigy[sigy], x0)
                        y_sim[x] = (data_x[x] + data_y[x]) - x_sim[x]
                    
                    # In case you need y(x): y_sim[x] = np.exp((grid_sigy[sigy]/grid_sigx[sigx])*np.log(x_sim[x]-x0) - grid_sigy[sigy]*np.log((grid_sigx[sigx] - 1)/(grid_sigy[sigy]-1)*(grid_sigy[sigy]/grid_sigx[sigx])))
                
                sqe[sigx,sigy,xi] = calculate_sqe(y_sim, data_x, data_y, x0_index, gridm, cons_bin, y_data_sel_sum, W)[0]
                
                # Not the first sigy evaluated
                if sigy > sigy_start:
                    # Checking if sigy passed minimum
                    if sqe[sigx,sigy,xi] > sqe[sigx,sigy-1,xi]:
                        sqe[sigx, :sigy_start, xi] = np.Inf
                        sqe[sigx, (sigy+1):, xi] = np.Inf
                        sigy_start = np.max(np.array([sigy - 2, np.int64(0)]))
                        break
                    # Checking if sigy did not reach minimum
                    elif sigy == len(grid_sigy):
                        sqe[sigx,:sigy_start, xi] = np.Inf
                        sigy_start = np.max(np.array([sigy - 1, np.int64(0)]))
    
                        
            progress = np.floor((sigx/len(grid_sigx))*100)
            if (progress / 20) > counter :
                counter = counter + 1
                print(progress, " % of x0 =",x0,"concluded")
    
    return sqe

# Calculating squared error of moments
@njit
def calculate_sqe(y_sim, data_x, data_y, x0_index, gridm, cons_bin, y_data_sel_sum, W):
    
    m_errors = np.zeros(shape = len(gridm))
    y_sim_sel_sum = np.zeros(shape = len(gridm))
    cons_bin_selected = cons_bin[:,x0_index]
    
    for b in range(len(m_errors)):
        y_sim_sel_sum[b] = np.sum(y_sim[cons_bin_selected == b])
    
    m_errors = ((y_sim_sel_sum - y_data_sel_sum[:,x0_index]) / y_data_sel_sum[:,x0_index])
    sqe = (m_errors.reshape((1,len(gridm))) @ W @ m_errors.reshape((len(gridm),1)))[0][0]
    
    return sqe, m_errors
        
        
@njit
def f(x, c, xi, sigma_x, sigma_y, x0):
    return x + (x-x0)**(sigma_y/sigma_x)*(((sigma_x-1)/sigma_x)*(sigma_y/(xi*(sigma_y-1))))**(-sigma_y) - c

@njit
def df(x, c, xi, sigma_x, sigma_y, x0):
    return 1 + sigma_y/sigma_x * (x-x0)**((sigma_y-sigma_x)/sigma_x)*(((sigma_x-1)/sigma_x)*(sigma_y/(xi*(sigma_y-1))))**(-sigma_y)

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

def create_cons_moment_grid(which):
    
    if which == "arbitrary":
        grid1 = np.arange(0, 1000*12, 200*12)
        grid2 = np.arange(1000*12, 2600*12, 400*12)
        grid3 = np.arange(2600*12, 10100*12,1500*12)
        grid4 = np.arange(10100*12, 20000*12, 5000*12)
        gridm = np.concatenate((grid1,grid2,grid3,grid4))
    
    if which == "percentile5":
        gridm = np.array([0,
                          144.36,
                          201.56,
                          251.27,
                          297.58,
                          345.78,
                          392.88,
                          445.24,
                          499.12,
                          556.54,
                          620.43,
                          688.46,
                          767.22,
                          859.76,
                          968.87,
                          1100.25,
                          1265.36,
                          1491.37,
                          1863.56,
                          2596.38])*12
    
    if which == "percentile10":
        gridm = np.array([0,
                          201.56,
                          297.58,
                          392.88,
                          499.12,
                          620.43,
                          767.22,
                          968.87,
                          1265.36,
                          1863.56,
                          2596.38])*12
        
    return gridm

@njit
def run_x0_simulations(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi, W):
    
    # Starting timer
    # start_time = time.time()
    print("Running initial calculations... \n")
    
    sqe = np.zeros(shape = (len(grid_x0), len(grid_sigx), len(grid_sigy), len(grid_xi)))
    consumpt, cons_bin, y_data_sel_sum, gridm = initial_computations(data_x, data_y, gridq, grid_x0)
    
    print("Proceeding to simulations... \n")
    for x0 in range(len(grid_x0)):
        
        sqe[x0,:,:,:] = run_simulations_by_x0(grid_sigx, grid_sigy, grid_x0[x0], grid_x0, data_x, data_y, gridm, grid_xi, consumpt, cons_bin, y_data_sel_sum, W)
        
        #hour = round(((time.time() - start_time)/60)//60)
        #minute = round((time.time() - start_time)//60 - hour*60)
        #second = round((time.time() - start_time) - hour*60*60 - minute*60)
        print("\n",np.floor((x0/len(grid_x0))*100), " % of all simulations concluded \n")             
        #print("Time elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")
    
    return sqe

@njit
def run_two_steps(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi):
    
    W = np.identity(len(gridq))
    
    print("Running first simulation...\n")
    
    sqe_I = run_x0_simulations(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi, W)
    
    print("Recalculating weight matrix...\n")
    
    sol = np.argwhere(sqe_I == np.min(sqe_I))
    sol_x0_1 = grid_x0[sol[0][0]]
    sol_sigx_1 = grid_sigx[sol[0][1]]
    sol_sigy_1 = grid_sigy[sol[0][2]]
    sol_xi_1 = grid_xi[sol[0][3]]
    errors_sol = simulate_single(sol_sigx_1, sol_sigy_1, sol_x0_1, data_x, data_y, gridq, sol_xi_1, W, grid_x0)[1]
    Omega, W_new = recalculate_weights(errors_sol)
    
    print(" First stage concluded!\n\n Parameters found: \n x0:",sol_x0_1,"\n sigx:",sol_sigx_1,"\n sigy:",sol_sigy_1,"\n xi:",sol_xi_1)
    
    print("Running second simulation...\n")
    
    sqe_W = run_x0_simulations(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi, W_new)
    
    sol = np.argwhere(sqe_W == np.min(sqe_W))
    sol_x0_2 = grid_x0[sol[0][0]]
    sol_sigx_2 = grid_sigx[sol[0][1]]
    sol_sigy_2 = grid_sigy[sol[0][2]]
    sol_xi_2 = grid_xi[sol[0][3]]
    
    print(" (First Stage Recap)\n\n Parameters found: \n x0:",sol_x0_1,"\n sigx:",sol_sigx_1,"\n sigy:",sol_sigy_1,"\n xi:",sol_xi_1,"\n---------------------------")
    print(" Second stage concluded!\n\n Parameters found: \n x0:",sol_x0_2,"\n sigx:",sol_sigx_2,"\n sigy:",sol_sigy_2,"\n xi:",sol_xi_2)

    return sqe_I, sqe_W, W_new, Omega

def select_data(pof_df):
    data_x = pof_df["non-temptation_pc"].to_numpy()
    data_y = pof_df["temptation_pc"].to_numpy()
    
    return data_x, data_y

@njit
def recalculate_weights(m_errors):
    Omega = (1/len(m_errors))* (m_errors.reshape((len(m_errors),1)) @ m_errors.reshape((1,len(m_errors))))
    W_new = np.linalg.pinv(Omega)
    
    return Omega, W_new

@njit
def simulate_single(sigx, sigy, x0, data_x, data_y, gridq, xi, W, grid_x0):
    
    consumpt, cons_bin, y_data_sel_sum, gridm = initial_computations(data_x, data_y, gridq, grid_x0)
    
    y_sim = np.zeros(shape = (len(data_x)))
    x_sim = np.zeros(shape = (len(data_x)))
    
    for x in range(len(data_x)):
        if consumpt[x] < (x0 + 2):
            x_sim[x] = consumpt[x]
        else:
            x_sim[x] = Newton(f, df, (x0 + 0.0001), 1e-6,int(1e6), consumpt[x], xi, sigx, sigy, x0)
            y_sim[x] = (data_x[x] + data_y[x]) - x_sim[x]
    
    x0_index = np.argwhere(grid_x0 == x0)[0][0]
    sqe, m_errors = calculate_sqe(y_sim, data_x, data_y, x0_index, gridm, cons_bin, y_data_sel_sum, W)
    
    return sqe, m_errors
    
@njit
def initial_computations(data_x, data_y, gridq, grid_x0):
    
    consumpt = data_x + data_y
    
    # Creating array to hold which consumption bin each consumption belongs to
    cons_bin = np.zeros(shape = (len(consumpt), len(grid_x0)))
    
    # Creating array with the sum of temptation consumption in each bin
    y_data_sel_sum = np.zeros(shape = (len(gridq), len(grid_x0)))
    
    # Creating grid of moments above x0
    gridm = np.zeros(shape = (len(gridq),len(grid_x0)))
    
    for x0 in range(len(grid_x0)):
        
        consumpt_censored = consumpt[consumpt>grid_x0[x0]]
        gridm[:,x0] = np.quantile(consumpt_censored,gridq)
        
        # Filling cons_bin accordingly
        for c in range(len(consumpt)): 
            cons_bin[c, x0] = np.int64(np.sum(gridm[:,x0] < consumpt[c]) - 1)
            
        for b in range(len(gridm)):
            y_data_sel_sum[b, x0] = np.sum(data_y[cons_bin[:, x0] == b])
        
    return consumpt, cons_bin, y_data_sel_sum, gridm