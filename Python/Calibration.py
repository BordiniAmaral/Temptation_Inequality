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
import time

def import_POF_data(file_path):
    pof_df = pd.read_csv(file_path, sep = ";")
    return pof_df

def select_data(pof_df):
    data_x = pof_df["non-temptation_pc"].to_numpy()
    data_y = pof_df["temptation_pc"].to_numpy()
    
    return data_x, data_y

@njit
def run_simulations_by_x0(grid_gam1, grid_gam2, x0, grid_x0, data_x, data_y, gridm, consumpt, cons_bin, y_data_sel_sum, W, verbose):
    
    if verbose: print("------------------------ \n   Current x0:",x0,"  \n------------------------\n")
    sqe = np.zeros(shape = (len(grid_gam1), len(grid_gam2)))
    counter = 0
    gam2_start = np.int32(0)
    
    for gam1 in range(len(grid_gam1)):
        for gam2 in range(gam2_start, len(grid_gam2)):
            
            y_sim = np.zeros(shape = (len(data_x)))
            x_sim = np.zeros(shape = (len(data_x)))
            
            for x in range(len(data_x)):
                if consumpt[x] < (x0 + 2):
                    x_sim[x] = consumpt[x]
                else:
                    x_sim[x] = Newton(f, df, (x0 + 1), 1e-6,int(1e6), consumpt[x], grid_gam1[gam1], grid_gam2[gam2], x0)
                    y_sim[x] = (data_x[x] + data_y[x]) - x_sim[x]
            
            sqe[gam1,gam2] = calculate_sqe(y_sim, data_x, data_y, gridm, cons_bin, y_data_sel_sum, W)[0]
        
            if gam2 > gam2_start+1:
                # Checking if gam2 passed minimum by 2 steps at least
                if sqe[gam1,gam2] > sqe[gam1,gam2-2]:
                    sqe[gam1,:gam2_start] = np.Inf
                    sqe[gam1,(gam2+1):] = np.Inf
                    gam2_start = np.max(np.array([gam2 - 3, np.int64(0)]))
                    break
                # Checking if sigy did not reach minimum
                elif gam2 == len(grid_gam2):
                    sqe[gam1,:gam2_start] = np.Inf
                    gam2_start = np.max(np.array([gam2 - 3, np.int64(0)]))
        
        progress = np.floor((gam1/len(grid_gam1))*100)
        if (progress / 20) > counter :
            counter = counter + 1
            if verbose: print(progress, " % of x0 =",x0,"concluded")
    
    return sqe

# Calculating squared error of moments
@njit
def calculate_sqe(y_sim, data_x, data_y, gridm, cons_bin, y_data_sel_sum, W):
    
    m_errors = np.zeros(shape = len(gridm))
    y_sim_sel_sum = np.zeros(shape = len(gridm))
    
    for b in range(len(m_errors)):
        y_sim_sel_sum[b] = np.sum(y_sim[cons_bin == b])
    
    m_errors = ((y_sim_sel_sum - y_data_sel_sum) / y_data_sel_sum)
    sqe = (m_errors.reshape((1,len(gridm))) @ W @ m_errors.reshape((len(gridm),1)))[0][0]
    
    return sqe, m_errors
        
        
@njit
def f(x, c, gam1, gam2, x0):
    return x + np.exp(gam1)*(x-x0)**gam2 - c

@njit
def df(x, c, gam1, gam2, x0):
    return 1 + np.exp(gam1)*gam2*(x-x0)**(gam2-1)

@njit
def Newton(f, df, x_init, e, maxit, c, gam1, gam2, x0):
    xn = x_init
    for n in range(0,maxit):
        fxn = f(xn, c, gam1, gam2, x0)
        if abs(fxn)<e:
            return xn
        dfxn = df(xn, c, gam1, gam2, x0)
        if dfxn == 0:
            print("Reached zero derivative in iteration",n,"please review initial guess")
            return None
        xn = xn - fxn/dfxn
    print("Maximum iterations (",n,") reached. \n Current x:",xn,"\n Current f(x):",fxn)
    print("\n C evaluated: ", c,"\n gam1 evaluated: ", gam1,"\n gam2 evaluated: ", gam2)
    return None

@njit
def run_x0_simulations(grid_gam1, grid_gam2, grid_x0, data_x, data_y, gridq, W, quantile, verbose):
    
    # Starting timer
    # start_time = time.time()
    if verbose: print("Running initial calculations... \n")
    
    sqe = np.zeros(shape = (len(grid_x0), len(grid_gam1), len(grid_gam2)))
    consumpt, cons_bin, y_data_sel_sum, gridm = initial_computations(data_x, data_y, gridq, quantile)
    
    if verbose: print("Proceeding to simulations... \n")
    for x0 in range(len(grid_x0)):
        
        sqe[x0,:,:] = run_simulations_by_x0(grid_gam1, grid_gam2, grid_x0[x0], grid_x0, data_x, data_y, gridm, consumpt, cons_bin, y_data_sel_sum, W, verbose)
        
        #hour = round(((time.time() - start_time)/60)//60)
        #minute = round((time.time() - start_time)//60 - hour*60)
        #second = round((time.time() - start_time) - hour*60*60 - minute*60)
        if verbose: print("\n",np.floor((x0/len(grid_x0))*100), " % of all simulations concluded \n")             
        #print("Time elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")
    
    return sqe
    
@njit
def initial_computations(data_x, data_y, gridq, quantile):
    
    consumpt = data_x + data_y
    
    # Creating array to hold which consumption bin each consumption belongs to
    cons_bin = np.zeros(shape = (len(consumpt)))
    
    # Creating array with the sum of temptation consumption in each bin
    y_data_sel_sum = np.zeros(shape = (len(gridq)))
    
    # Creating grid of moments
    if quantile:
        gridm = np.quantile(consumpt,gridq)
    if not quantile:
        gridm = gridq
    
    for c in range(len(consumpt)): 
        cons_bin[c] = np.int64(np.sum(gridm <= consumpt[c]) - 1)
            
    for b in range(len(gridm)):
        y_data_sel_sum[b] = np.sum(data_y[cons_bin == b])
    
        
    return consumpt, cons_bin, y_data_sel_sum, gridm

def bootstrap(samples, sample_size, grid_x0, grid_gam1, grid_gam2, data_x, data_y, gridq, W, quantile, verbose):
    
    start_time = time.time()
    
    boot_estimates = np.zeros(shape = (samples,3))
    corner_solution = np.zeros(shape = (samples,3))
    pop_size = len(data_x)
    
    print("\nRunning bootstrap... this will take some time, relax")
    
    for s in range(samples):
        
        indexes_sampled = np.random.choice(pop_size, size = sample_size)
        data_x_sampled = data_x[indexes_sampled]
        data_y_sampled = data_y[indexes_sampled]
        
        sqe = run_x0_simulations(grid_gam1, grid_gam2, grid_x0, data_x_sampled, data_y_sampled, gridq, W, quantile, verbose)
        sol_index = np.argwhere(sqe == np.min(sqe))
        boot_estimates[s,:] = np.array([grid_x0[sol_index[0][0]],grid_gam1[sol_index[0][1]],grid_gam2[sol_index[0][2]]])
        corner_solution[s,:] = np.array([int(sol_index[0][0]==(len(grid_x0)-1))-int(sol_index[0][0]==0), \
                                         int(sol_index[0][1]==(len(grid_gam1)-1))-int(sol_index[0][1]==0), \
                                         int(sol_index[0][2]==(len(grid_gam2)-1))-int(sol_index[0][2]==0)])
        
        print(s+1,"out of",samples,"calculated")
    
    IC95_x0 = np.quantile(boot_estimates[:,0],[0.025,0.975])
    IC95_gam1 = np.quantile(boot_estimates[:,1],[0.025,0.975])
    IC95_gam2 = np.quantile(boot_estimates[:,2],[0.025,0.975])
    
    avg_x0 = np.sum(boot_estimates[:,0])/samples
    avg_gam1 = np.sum(boot_estimates[:,1])/samples
    avg_gam2 = np.sum(boot_estimates[:,2])/samples
    
    viol_x0 = np.array([np.sum(corner_solution[:,0]==-1),np.sum(corner_solution[:,0]==1)])
    viol_gam1 = np.array([np.sum(corner_solution[:,1]==-1),np.sum(corner_solution[:,1]==1)])
    viol_gam2 = np.array([np.sum(corner_solution[:,2]==-1),np.sum(corner_solution[:,2]==1)])
    
    end_time = time.time()
    hour = round(((end_time - start_time)/60)//60)
    minute = round((end_time - start_time)//60 - hour*60)
    second = round((end_time - start_time) - hour*60*60 - minute*60)
    print("Bootstrap finished! \nTime elapsed (Total): ", hour, " h ", minute, "min ", second, "s \n \n")
    
    print("----------------------------------------------")
    print("             Reporting values")
    print("----------------------------------------------")
    print("x0: \n  Mean = ",np.round(avg_x0,4)/12,"\n  IC95 = [",IC95_x0[0]/12,",",IC95_x0[1]/12,"]","\n  Space = [",grid_x0[0]/12,",",grid_x0[-1]/12,"]","\n  Violations:",viol_x0[0],"lower,",viol_x0[1],"upper")
    print("----------------------------------------------")
    print("gam1: \n  Mean = ",np.round(avg_gam1,4),"\n  IC95 = [",np.round(IC95_gam1[0],4),",",np.round(IC95_gam1[1],4),"]","\n  Space = [",grid_gam1[0],",",grid_gam1[-1],"]","\n  Violations:",viol_gam1[0],"lower,",viol_gam1[1],"upper")
    print("----------------------------------------------")
    print("gam2: \n  Mean = ",np.round(avg_gam2,4),"\n  IC95 = [",np.round(IC95_gam2[0],4),",",np.round(IC95_gam2[1],4),"]","\n  Space = [",grid_gam2[0],",",grid_gam2[-1],"]","\n  Violations:",viol_gam2[0],"lower,",viol_gam2[1],"upper")
    print("----------------------------------------------")
    
    return boot_estimates, corner_solution, IC95_x0, IC95_gam1, IC95_gam2

@njit
def compute_sim_tempt_frac(data_x, data_y, x0, gam1, gam2, gridq, quantile, verbose):
    
    y_sim = np.zeros(shape = (len(data_x)))
    x_sim = np.zeros(shape = (len(data_x)))
    
    consumpt, cons_bin, y_data_sel_sum, gridm = initial_computations(data_x, data_y, gridq, quantile)
    
    for x in range(len(data_x)):
        if consumpt[x] < (x0 + 2):
            x_sim[x] = consumpt[x]
        else:
            x_sim[x] = Newton(f, df, (x0 + 1), 1e-6,int(1e6), consumpt[x], gam1, gam2, x0)
            y_sim[x] = consumpt[x] - x_sim[x]
            
    fraction_sim = np.zeros(shape = len(gridm))
    
    for b in range(len(gridm)):
        y_sim_sel_sum = np.sum(y_sim[cons_bin == b])
        cons_sel_sum = np.sum(consumpt[cons_bin == b])
        
        fraction_sim[b] = y_sim_sel_sum / cons_sel_sum
        
    if verbose:
        print("----------------------------------------------")
        print("             Reporting values")
        print("----------------------------------------------")
        print("     Region          Sim.Tempt. Fraction ")
        for b in range(len(gridm)):
            print("    ",b,"              ",np.round(fraction_sim[b],4)*100,"%")
    
    return fraction_sim
        