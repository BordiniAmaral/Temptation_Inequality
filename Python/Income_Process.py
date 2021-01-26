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