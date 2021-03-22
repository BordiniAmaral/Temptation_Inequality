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
#                             Graphics                                       #
##############################################################################

# This file is dedicated to generating the desired plots

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import make_interp_spline

def capital_distributions(n, grida, gridz, distr_mass, k_mass):
    
    lifecycle_mass = np.sum(distr_mass[:n,:,:])
    mass_by_k = np.zeros(shape = len(grida))
    mass_by_age_k = np.zeros(shape = (n, len(grida)))
    
    for a in range(len(grida)):
        mass_by_k[a] = np.sum(distr_mass[:n,:,a]) / lifecycle_mass
        for period in range(n):
            mass_by_age_k[period,a] = np.sum(distr_mass[period,:,a]) / lifecycle_mass
    
    return mass_by_k, mass_by_age_k

def plot_total_k_distr(mass_by_k, grida, bin_size, description):
    
    # Detecting first and last positive mass
    first_index = np.argwhere(mass_by_k > 0)[0][0]
    last_index = np.argwhere(mass_by_k > 0)[-1][0]
    
    if grida[first_index] < 0:
        first_value = (grida[first_index] // bin_size ) * bin_size
    else:
        first_value = grida[first_index]
    
    a_plot_grid = np.arange(np.min([0, first_value]),grida[last_index], bin_size)
    a_plot_values = np.zeros(shape = len(a_plot_grid))
    
    for a in range(len(a_plot_grid)):
        a_plot_values[a] = np.sum(mass_by_k[abs(grida - a_plot_grid[a]) < (bin_size/2)])
    
    spl = make_interp_spline(a_plot_grid[a_plot_values>0],a_plot_values[a_plot_values>0], k = 3)
    smooth_a = spl(a_plot_grid)
    # Plotting graph for total k distribution
    
    plt.figure(figsize=(7,7))
    plt.fill_between(a_plot_grid, 0, smooth_a)
    plt.ylabel('Mass of agents')
    plt.xlabel('Asset level')
    plt.title('Capital Distribution - '+ description)

def plot_total_k_lorenz(mass_by_k, grida, description):
    
    cum_mass_by_k = np.cumsum(mass_by_k)
    cum_kfrac_by_k = np.cumsum(mass_by_k * grida) / np.sum(mass_by_k * grida)
    
    plt.figure(figsize=(7,7))
    plt.fill_between(cum_mass_by_k, 0, cum_kfrac_by_k, color = 'lightsteelblue')
    plt.plot(cum_mass_by_k, cum_mass_by_k, linestyle = 'dashed', color='black')
    plt.ylabel('Cumulative Wealth')
    plt.xlabel('Cumulative Population')
    plt.title('Lorenz Curve - '+ description)
    
def compare_total_k_distr(mass_by_k1, mass_by_k2, mass_by_k3, grida, bin_size, description, label1, label2, label3, log = False, trim_upper = False, trim_value = 0):
    
    # Detecting first and last positive mass
    first_index1 = np.argwhere(mass_by_k1 > 0)[0][0]
    last_index1 = np.argwhere(mass_by_k1 > 0)[-1][0]
    
    first_index2 = np.argwhere(mass_by_k2 > 0)[0][0]
    last_index2 = np.argwhere(mass_by_k2 > 0)[-1][0]
    
    first_index3 = np.argwhere(mass_by_k3 > 0)[0][0]
    last_index3 = np.argwhere(mass_by_k3 > 0)[-1][0]
    
    if (grida[first_index1] < 0 or grida[first_index2] < 0 or grida[first_index3] < 0) and not log:
        first_value = (grida[min([first_index1, first_index2,first_index3])] // bin_size ) * bin_size
        last_value = grida[max([last_index1,last_index2,last_index3])]
        if trim_upper:
            last_value = min(last_value, trim_value)
    elif not log:
        first_value = grida[min([first_index1, first_index2,first_index3])]
        last_value = grida[max([last_index1,last_index2,last_index3])]
        if trim_upper:
            last_value = min(last_value, trim_value)
    else:
        first_value = np.log(grida[np.argwhere(grida == 0)[0][0]+1])
        last_value = np.log(grida[max([last_index1,last_index2,last_index3])])
    
    
    if log:
        mass_by_k1 = mass_by_k1[np.argwhere(grida == 0)[0][0]+1:]
        mass_by_k2 = mass_by_k2[np.argwhere(grida == 0)[0][0]+1:]
        mass_by_k3 = mass_by_k2[np.argwhere(grida == 0)[0][0]+1:]
        
        first_index1 = np.argwhere(grida == 0)[0][0]+1
        first_index2 = np.argwhere(grida == 0)[0][0]+1
        first_index3 = np.argwhere(grida == 0)[0][0]+1
        
        a_plot_grid = np.arange(first_value, last_value+bin_size, bin_size)
        
        a_plot_values1 = np.zeros(shape = len(a_plot_grid))
        a_plot_values2 = np.zeros(shape = len(a_plot_grid))
        a_plot_values3 = np.zeros(shape = len(a_plot_grid))
        
        for a in range(len(a_plot_grid)):
            a_plot_values1[a] = np.sum(mass_by_k1[abs(np.log(grida[grida>0]) - a_plot_grid[a]) < (bin_size/2)])
            a_plot_values2[a] = np.sum(mass_by_k2[abs(np.log(grida[grida>0]) - a_plot_grid[a]) < (bin_size/2)])
            a_plot_values3[a] = np.sum(mass_by_k3[abs(np.log(grida[grida>0]) - a_plot_grid[a]) < (bin_size/2)])
        
    else:
        a_plot_grid = np.arange(first_value, last_value, bin_size)
        
        a_plot_values1 = np.zeros(shape = len(a_plot_grid))
        a_plot_values2 = np.zeros(shape = len(a_plot_grid))
        a_plot_values3 = np.zeros(shape = len(a_plot_grid))
        
        for a in range(len(a_plot_grid)):
            a_plot_values1[a] = np.sum(mass_by_k1[abs(grida - a_plot_grid[a]) < (bin_size/2)])
            a_plot_values2[a] = np.sum(mass_by_k2[abs(grida - a_plot_grid[a]) < (bin_size/2)])
            a_plot_values3[a] = np.sum(mass_by_k3[abs(grida - a_plot_grid[a]) < (bin_size/2)])
            
            
    spl1 = make_interp_spline(a_plot_grid[a_plot_values1>0],a_plot_values1[a_plot_values1>0], k = 3)
    spl2 = make_interp_spline(a_plot_grid[a_plot_values2>0],a_plot_values2[a_plot_values2>0], k = 3)
    spl3 = make_interp_spline(a_plot_grid[a_plot_values3>0],a_plot_values2[a_plot_values3>0], k = 3)
    smooth_a1 = spl1(a_plot_grid)
    smooth_a2 = spl2(a_plot_grid)
    smooth_a3 = spl3(a_plot_grid)
    
    # Plotting graph for total k distribution
    
    plt.figure(figsize=(7,7))
    plt.fill_between(a_plot_grid, 0, smooth_a1, alpha = 0.4, label = label1)
    plt.fill_between(a_plot_grid, 0, smooth_a2, alpha = 0.4, label = label2)
    plt.fill_between(a_plot_grid, 0, smooth_a3, alpha = 0.4, label = label3)
    plt.ylabel('Mass of agents')
    plt.xlabel('Log '*log + 'Asset level')
    plt.title('Capital Distribution - '+ description)
    plt.legend(title = "Distributions")
    
def compare_total_k_lorenz(mass_by_k1, mass_by_k2, mass_by_k3, grida, description, label1, label2, label3):
    
    cum_mass_by_k1 = np.cumsum(mass_by_k1)
    cum_kfrac_by_k1 = np.cumsum(mass_by_k1 * grida) / np.sum(mass_by_k1 * grida)
    
    cum_mass_by_k2 = np.cumsum(mass_by_k2)
    cum_kfrac_by_k2 = np.cumsum(mass_by_k2 * grida) / np.sum(mass_by_k2 * grida)
    
    cum_mass_by_k3 = np.cumsum(mass_by_k3)
    cum_kfrac_by_k3 = np.cumsum(mass_by_k3 * grida) / np.sum(mass_by_k3 * grida)
    
    plt.figure(figsize=(7,7))
    plt.plot(cum_mass_by_k1, cum_kfrac_by_k1, label = label1)
    plt.plot(cum_mass_by_k2, cum_kfrac_by_k2, label = label2)
    plt.plot(cum_mass_by_k3, cum_kfrac_by_k3, label = label3)
    plt.plot(cum_mass_by_k1, cum_mass_by_k1, linestyle = 'dashed', color='black')
    plt.ylabel('Cumulative Wealth')
    plt.xlabel('Cumulative Population')
    plt.title('Lorenz Curve - '+ description)
    plt.legend(title = "Distributions")
    
def plot_k_evolution(age_start,n,mass_by_age_k, quants, grida, description):
    
    cumulative_per_age = np.cumsum(mass_by_age_k,axis = 1) * (n+1)
    
    # Finding indexes for each quantile at each age
    quant_index = np.zeros(shape = (len(quants),n))
    quant_value = np.zeros(shape = (len(quants),n))
    for age in range(n):
        for q in range(len(quants)):
            quant_index[q,age] = np.sum(cumulative_per_age[age,:] < quants[q])
            quant_value[q,age] = grida[int(quant_index[q,age])-1]
    
    # Plotting
    age_tick = np.arange(age_start,age_start+n)
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    for q in range(len(quants)):
        ax.plot(age_tick,quant_value[q,:], label = "Quantile "+str(quants[q]))
    ax.legend(loc='upper left')
    fig.suptitle('Asset evolution by quantile - '+description)
    ax.set_xlabel('Household head age')
    ax.set_ylabel('Household asset position [BRL]')
    ax.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.show()
    
def compare_k_evolution(age_start, n, mass_by_age_k1, mass_by_age_k2, mass_by_age_k3, quants, grida, description1, description2, description3):
    
    cohort_mass = np.sum(mass_by_age_k1[0,:])
    cumulative_per_age1 = np.round(np.cumsum(mass_by_age_k1,axis = 1) / cohort_mass , 4)
    cumulative_per_age2 = np.round(np.cumsum(mass_by_age_k2,axis = 1) / cohort_mass , 4)
    cumulative_per_age3 = np.round(np.cumsum(mass_by_age_k3,axis = 1) / cohort_mass , 4)
    
    # Finding indexes for each quantile at each age
    quant_index1 = np.zeros(shape = (len(quants),n))
    quant_index2 = np.zeros(shape = (len(quants),n))
    quant_index3 = np.zeros(shape = (len(quants),n))
    quant_value1 = np.zeros(shape = (len(quants),n))
    quant_value2 = np.zeros(shape = (len(quants),n))
    quant_value3 = np.zeros(shape = (len(quants),n))
    for age in range(n):
        for q in range(len(quants)):
            quant_index1[q,age] = np.sum(cumulative_per_age1[age,:] < quants[q])
            quant_index2[q,age] = np.sum(cumulative_per_age2[age,:] < quants[q])
            quant_index3[q,age] = np.sum(cumulative_per_age3[age,:] < quants[q])
            quant_value1[q,age] = grida[int(quant_index1[q,age])-1]
            quant_value2[q,age] = grida[int(quant_index2[q,age])-1]
            quant_value3[q,age] = grida[int(quant_index3[q,age])-1]
    
    # Plotting
    age_tick = np.arange(age_start,age_start+n)
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(quants))))
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    for q in range(len(quants)):
        c = next(color)
        label = str(np.int(quants[q]*100)) + "%"
        ax.plot(age_tick,quant_value1[q,:], color = c, linestyle = '-', label = label)
        ax.plot(age_tick,quant_value2[q,:], color = c, linestyle = '--')
        ax.plot(age_tick,quant_value3[q,:], color = c, linestyle = '-.')
        
    handles, labels = ax.get_legend_handles_labels()
    display = list(range(0,len(quants)))
    Artist1 = plt.Line2D((0,1),(0,0), color='k', linestyle='-')
    Artist2 = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
    Artist3 = plt.Line2D((0,1),(0,0), color='k', linestyle='-.')
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[Artist1,Artist2,Artist3], \
              [label for i,label in enumerate(labels) if i in display]+[description1, description2, description3], loc='upper left')
    
    fig.suptitle('Wealth evolution by wealth quantile (per capita)')
    ax.set_xlabel('Household head age')
    ax.set_ylabel('Household wealth per capita [BRL]')
    ax.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.show()

def savings_and_income(n, grida, choice_a, gridz, r, include_interest):
    
    # Amount saved/dissaved by each state (age, productivity, asset position)
    # Considering savings as change in asset position
    
    gross_savings = np.zeros(shape = (n,len(gridz[:,0]),len(grida)))
    savings_rate = np.zeros(shape = (n,len(gridz[:,0]),len(grida)))
    total_income = np.zeros(shape = (n,len(gridz[:,0]),len(grida)))
    
    for age in range(n):
        for z in range(len(gridz[:,0])):
            for a in range(len(grida)):
                gross_savings[age,z,a] = grida[int(choice_a[age,z,a])] - grida[a]
                total_income[age,z,a] = gridz[z,age] + r*grida[a]*include_interest
                savings_rate[age,z,a] = gross_savings[age,z,a] / total_income[age,z,a]
    
    return gross_savings, savings_rate, total_income

def savings_by_quants(n, grida, choice_a, gridz, r, distr_mass, quants, include_interest):
    
    gross_savings, savings_rate, total_income = savings_and_income(n, grida, choice_a, gridz, r, include_interest)
    
    lifecycle_mass = np.sum(distr_mass[:n,:,:])
    
    quant_index = np.zeros(shape = (len(quants)))
    quant_mean = np.zeros(shape = (len(quants)))    # average savings rate by quant
    quant_wt_sd = np.zeros(shape = (len(quants)))   # weighted standard deviation of savings by quant
    
    # Turning multi-dimensional arrays into 1d
    flat_income = np.ndarray.flatten(total_income[:n,:,:])
    flat_mass = np.ndarray.flatten(distr_mass[:n,:,:])
    flat_savings = np.ndarray.flatten(gross_savings[:n,:,:])
    
    # Finding the increasing order of total income
    flat_income_order = np.argsort(flat_income)
    ordered_income = np.sort(flat_income)
    
    # Ordering income and respective masses
    mass_ordered_by_income = flat_mass[flat_income_order]
    gross_savings_ordered_by_income = flat_savings[flat_income_order]
    
    # Getting the cumulative mass in each age, ordered by income
    cumulative_mass = np.cumsum(mass_ordered_by_income) / lifecycle_mass
    
    for q in range(len(quants)):
        quant_index[q] = np.sum(cumulative_mass < quants[q])
        
        if q == 0:
            low = int(0)
        else:
            low = int(quant_index[q-1])
        top = int(quant_index[q])
            
        quant_mean[q] = np.sum(mass_ordered_by_income[low:top] * gross_savings_ordered_by_income[low:top]) \
            / np.sum(mass_ordered_by_income[low:top] * ordered_income[low:top])
        quant_wt_sd[q] = np.sqrt(np.sum((ordered_income[low:top]*mass_ordered_by_income[low:top])/np.sum(mass_ordered_by_income[low:top] * ordered_income[low:top]) \
            * (gross_savings_ordered_by_income[low:top]/ordered_income[low:top] - quant_mean[q])**2))
    
    return quant_mean, quant_wt_sd

def savings_rate_by_quants_age(n, grida, choice_a, gridz, r, distr_mass, quants, include_interest):
    
    gross_savings, savings_rate, total_income = savings_and_income(n, grida, choice_a, gridz, r, include_interest)
    
    cohort_mass = np.sum(distr_mass[0,:,:])
    
    quant_age_index = np.zeros(shape = (len(quants),n))
    quant_age_value = np.zeros(shape = (len(quants),n)) # average savings rate by quant and age
    
    for age in range(n):
        
        # Turning multi-dimensional arrays into 1d
        flat_income = np.ndarray.flatten(total_income[age,:,:])
        flat_mass = np.ndarray.flatten(distr_mass[age,:,:])
        flat_savings = np.ndarray.flatten(gross_savings[age,:,:])
        
        # Finding the increasing order of total income
        flat_income_order = np.argsort(flat_income)
        ordered_income = np.sort(flat_income)
        
        # Ordering income and respective masses
        mass_ordered_by_income = flat_mass[flat_income_order]
        gross_savings_ordered_by_income = flat_savings[flat_income_order]
        
        # Getting the cumulative mass in each age, ordered by income
        cumulative_mass = np.round(np.cumsum(mass_ordered_by_income) / cohort_mass, 5)
        
        for q in range(len(quants)):
            quant_age_index[q,age] = np.sum(cumulative_mass < quants[q])
            
            if q == 0:
                low = int(0)
            else:
                low = int(quant_age_index[q-1,age])
            top = int(quant_age_index[q,age])
                
            quant_age_value[q,age] = np.sum(mass_ordered_by_income[low:top] * gross_savings_ordered_by_income[low:top]) \
                / np.sum(mass_ordered_by_income[low:top] * ordered_income[low:top])
    
    return quant_age_value

def plot_savings_rate(age_start, n, gross_savings, savings_rate, total_income, distr_mass, quants, description, include_interest):
    
    quant_value = savings_rate_by_quants_age(n, gross_savings, total_income, distr_mass, quants, include_interest)
    
     # Plotting
    age_tick = np.arange(age_start,age_start+n)
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    for q in range(len(quants)):
        if q == 0:
            low = str(0)
        else:
            low = str(np.int(quants[q-1]*100))
        top = str(np.int(quants[q]*100))
        label = low + " - " + top + "%"
        ax.plot(age_tick, np.round(quant_value[q,:], decimals = 2), label = label)
    ax.plot(age_tick,np.repeat(0,n), linestyle = '--', color= 'black', linewidth=0.9)
    ax.legend(loc='lower left')
    fig.suptitle('Savings rate by total income quantile - '+description)
    ax.set_xlabel('Household head age')
    ax.set_ylabel('Average savings rate')
    fig.show()

def compare_savings_rate(age_start, n, n_select, quants, grida, gridz, r1, r2, r3, choice_a1, choice_a2, choice_a3, distr_mass1, distr_mass2, distr_mass3, include_interest, description1, description2, description3):
    
    quant_value1 = savings_rate_by_quants_age(n, grida, choice_a1, gridz, r1, distr_mass1, quants, include_interest)
    quant_value2 = savings_rate_by_quants_age(n, grida, choice_a2, gridz, r2, distr_mass2, quants, include_interest)
    quant_value3 = savings_rate_by_quants_age(n, grida, choice_a3, gridz, r3, distr_mass3, quants, include_interest)
    
     # Plotting
    age_tick = np.arange(age_start,age_start+n_select)
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(quants))))
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    for q in range(len(quants)):
        c = next(color)
        if q == 0:
            low = str(0)
        else:
            low = str(np.int(quants[q-1]*100))
        top = str(np.int(quants[q]*100))
        label = low + " - " + top + "%"
        ax.plot(age_tick, np.round(quant_value1[q,:n_select], decimals = 2), label = label, color = c, linestyle = "-" )
        ax.plot(age_tick, np.round(quant_value2[q,:n_select], decimals = 2), color = c, linestyle = "--")
        ax.plot(age_tick, np.round(quant_value3[q,:n_select], decimals = 2), color = c, linestyle = "-.")
    ax.plot(age_tick,np.repeat(0,n_select), linestyle = '--', color= 'black', linewidth=0.9)
    
    handles, labels = ax.get_legend_handles_labels()
    display = list(range(1,len(quants)))
    Artist1 = plt.Line2D((0,1),(0,0), color='k', linestyle='-')
    Artist2 = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
    Artist3 = plt.Line2D((0,1),(0,0), color='k', linestyle='-.')
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[Artist1,Artist2,Artist3], \
              [label for i,label in enumerate(labels) if i in display]+[description1, description2, description3], loc='lower left')
    
    fig.suptitle('Savings rate by total income quantile')
    ax.set_xlabel('Household head age')
    ax.set_ylabel('Average savings rate')
    fig.show()
    
    return quant_value1, quant_value2

def savings_and_wealth_report(n, mass_by_k, grida, quants, distr_mass, show_zero):
    
    zero_asset_index = np.int64(np.argwhere(grida == 0)[0][0])
    mass_zero = mass_by_k[zero_asset_index]
    
    cum_frac_by_k = np.cumsum(mass_by_k) / np.sum(mass_by_k)
    cum_kfrac_by_k = np.cumsum(mass_by_k * grida) / np.sum(mass_by_k * grida)
    
    # Finding indexes for each quantile
    quant_mass = np.zeros(shape = (len(quants)))
    
    for q in range(len(quants)):
        quant_index = np.sum(cum_frac_by_k < quants[q])
        quant_mass[q] = 1 - cum_kfrac_by_k[int(quant_index)]
    
    # Reporting Stats
    print("       ---------------------------\n          Wealth Distribution   \n       ---------------------------")
    if show_zero:
        print("\nHouseholds with zero or negative asset: ", np.round(mass_zero*100,1),"%")
    print("\n Richest        Share of Wealth")
    for q in range(len(quants)):
        print(" ",np.round((1 - quants[q])*100,0),"%          ",np.round(quant_mass[q]*100,1))
    
    
def plot_k_curves(k_supply1 ,k_supply2, k_supply3, grid_r, k_demand, label1, label2, label3):
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(grid_r,k_demand, label = "Demand curve", linestyle = '-', color = 'black')
    ax.plot(grid_r,k_supply1, label = label1, linestyle = '--')
    ax.plot(grid_r,k_supply2, label = label2, linestyle = '--')
    ax.plot(grid_r,k_supply3, label = label3, linestyle = '--')
    ax.legend(loc='upper right')
    fig.suptitle('Capital demand and supply curves')
    ax.set_xlabel('Net interest rate')
    ax.set_ylabel('Capital')
    plt.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    fig.show()