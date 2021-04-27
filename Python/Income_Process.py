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
#                          Income Process                                    #
##############################################################################

# Importing income process
import numpy as np

def load_income_process(import_path, files):
    
    import pandas as pd
    
    import_files = [import_path + "\\" + i for i in files]

    transition = pd.read_csv(import_files[0], sep = ",")
    values = pd.read_csv(import_files[2], sep = ";")
    
    # Adapting values into an array
    
    values_reshape = values.pivot(index = 'quant_v', columns = 'hh_age_head', values = 'value')
    values_array = values_reshape.to_numpy()
    
    transition_array = transition.to_numpy()
    
    return(transition_array, values_array)

def load_test_transition(mass_z, persist, width):
    
    Pi = np.zeros(shape = (len(mass_z),len(mass_z)))
    step = np.zeros(shape = int((width-1)/2))
    
    unif = 1/width
    mid = unif + (1 - unif)*persist
    
    for i in range(len(step)):
        if i == 0:
            step[i] = 1
        else:
            step[i] = 2*step[i-1]
    
    step = step / np.sum(step) 
    step = step * (1 - mid) / 2
    spread = np.concatenate((step, [mid] ,np.flip(step)))
    
    for i in range(len(mass_z)):
        
        if i < (width-1)/2:
            Pi[i,0] = np.sum(spread[0:int((width-1)/2+1-i)])
            Pi[i,1:int(width-(width-1)/2+i)] = spread[int((width-1)/2+1-i):]
        elif i >= (len(mass_z) - (width-1)/2):
            Pi[i,-1] = np.sum(spread[int(len(mass_z)- i + (width-1)/2-1):])
            Pi[i,int(i-(width-1)/2):-1] = spread[:int(len(mass_z)-i+(width-1)/2-1)]
        else:
            Pi[i,int(i - (width-1)/2):int(i+(width-1)/2+1)] = spread[:]
    
    return Pi