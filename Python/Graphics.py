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
    
    mass_by_k = np.zeros(shape = len(grida))
    mass_by_age_k = np.zeros(shape = (n, len(grida)))
    
    for a in range(len(grida)):
        mass_by_k[a] = np.sum(distr_mass[:,:,a])
        for period in range(n):
            mass_by_age_k[period,a] = np.sum(distr_mass[period,:,a])
    
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
    
def compare_total_k_distr(mass_by_k1, mass_by_k2, grida, bin_size, description, label1, label2, log = False, trim_upper = False, trim_value = 0):
    
    # Detecting first and last positive mass
    first_index1 = np.argwhere(mass_by_k1 > 0)[0][0]
    last_index1 = np.argwhere(mass_by_k1 > 0)[-1][0]
    
    first_index2 = np.argwhere(mass_by_k2 > 0)[0][0]
    last_index2 = np.argwhere(mass_by_k2 > 0)[-1][0]
    
    if (grida[first_index1] < 0 or grida[first_index2] < 0) and not log:
        first_value = (grida[min([first_index1, first_index2])] // bin_size ) * bin_size
        last_value = grida[max([last_index1,last_index2])]
        if trim_upper:
            last_value = min(last_value, trim_value)
    elif not log:
        first_value = grida[min([first_index1, first_index2])]
        last_value = grida[max([last_index1,last_index2])]
        if trim_upper:
            last_value = min(last_value, trim_value)
    else:
        first_value = np.log(grida[np.argwhere(grida == 0)[0][0]+1])
        last_value = np.log(grida[max([last_index1,last_index2])])
    
    
    if log:
        mass_by_k1 = mass_by_k1[np.argwhere(grida == 0)[0][0]+1:]
        mass_by_k2 = mass_by_k2[np.argwhere(grida == 0)[0][0]+1:]
        
        first_index1 = np.argwhere(grida == 0)[0][0]+1
        first_index2 = np.argwhere(grida == 0)[0][0]+1
        
        a_plot_grid = np.arange(first_value, last_value+bin_size, bin_size)
        
        a_plot_values1 = np.zeros(shape = len(a_plot_grid))
        a_plot_values2 = np.zeros(shape = len(a_plot_grid))
        
        for a in range(len(a_plot_grid)):
            a_plot_values1[a] = np.sum(mass_by_k1[abs(np.log(grida[grida>0]) - a_plot_grid[a]) < (bin_size/2)])
            a_plot_values2[a] = np.sum(mass_by_k2[abs(np.log(grida[grida>0]) - a_plot_grid[a]) < (bin_size/2)])
        
    else:
        a_plot_grid = np.arange(first_value, last_value, bin_size)
        
        a_plot_values1 = np.zeros(shape = len(a_plot_grid))
        a_plot_values2 = np.zeros(shape = len(a_plot_grid))
        
        for a in range(len(a_plot_grid)):
            a_plot_values1[a] = np.sum(mass_by_k1[abs(grida - a_plot_grid[a]) < (bin_size/2)])
            a_plot_values2[a] = np.sum(mass_by_k2[abs(grida - a_plot_grid[a]) < (bin_size/2)])
            
    spl1 = make_interp_spline(a_plot_grid[a_plot_values1>0],a_plot_values1[a_plot_values1>0], k = 3)
    spl2 = make_interp_spline(a_plot_grid[a_plot_values2>0],a_plot_values2[a_plot_values2>0], k = 3)
    smooth_a1 = spl1(a_plot_grid)
    smooth_a2 = spl2(a_plot_grid)
    
    # Plotting graph for total k distribution
    
    plt.figure(figsize=(7,7))
    plt.fill_between(a_plot_grid, 0, smooth_a1, alpha = 0.4, label = label1)
    plt.fill_between(a_plot_grid, 0, smooth_a2, alpha = 0.4, label = label2)
    plt.ylabel('Mass of agents')
    plt.xlabel('Log '*log + 'Asset level')
    plt.title('Capital Distribution - '+ description)
    plt.legend(title = "Distributions")
    
def compare_total_k_lorenz(mass_by_k1, mass_by_k2, grida, description, label1, label2):
    
    cum_mass_by_k1 = np.cumsum(mass_by_k1)
    cum_kfrac_by_k1 = np.cumsum(mass_by_k1 * grida) / np.sum(mass_by_k1 * grida)
    
    cum_mass_by_k2 = np.cumsum(mass_by_k2)
    cum_kfrac_by_k2 = np.cumsum(mass_by_k2 * grida) / np.sum(mass_by_k2 * grida)
    
    plt.figure(figsize=(7,7))
    plt.plot(cum_mass_by_k1, cum_kfrac_by_k1, label = label1)
    plt.plot(cum_mass_by_k2, cum_kfrac_by_k2, label = label2)
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
    
def compare_k_evolution(age_start,n,mass_by_age_k1, mass_by_age_k2, quants, grida, description1, description2):
    
    cumulative_per_age1 = np.cumsum(mass_by_age_k1,axis = 1) * n
    cumulative_per_age2 = np.cumsum(mass_by_age_k2,axis = 1) * n
    
    # Finding indexes for each quantile at each age
    quant_index1 = np.zeros(shape = (len(quants),n))
    quant_index2 = np.zeros(shape = (len(quants),n))
    quant_value1 = np.zeros(shape = (len(quants),n))
    quant_value2 = np.zeros(shape = (len(quants),n))
    for age in range(n):
        for q in range(len(quants)):
            quant_index1[q,age] = np.sum(cumulative_per_age1[age,:] < quants[q])
            quant_index2[q,age] = np.sum(cumulative_per_age2[age,:] < quants[q])
            quant_value1[q,age] = grida[int(quant_index1[q,age])-1]
            quant_value2[q,age] = grida[int(quant_index2[q,age])-1]
    
    # Plotting
    age_tick = np.arange(age_start,age_start+n)
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    for q in range(len(quants)):
        ax.plot(age_tick,quant_value1[q,:], label = description1+" "+str(quants[q]),linestyle = '--')
        ax.plot(age_tick,quant_value2[q,:], label = description2+" "+str(quants[q]))
    ax.legend(loc='upper left')
    fig.suptitle('Asset evolution by quantile')
    ax.set_xlabel('Household head age')
    ax.set_ylabel('Household asset position [BRL]')
    ax.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.show()

def savings_rate(n, grida, choice_a, gridz, mass_z, r):
    
    # Amount saved/dissaved by each state (age, productivity, asset position)
    # Considering savings as change in asset position
    
    gross_savings = np.zeros(shape = (n,len(mass_z),len(grida)))
    savings_rate = np.zeros(shape = (n,len(mass_z),len(grida)))
    
    for age in range(n):
        for z in range(len(mass_z)):
            for a in range(len(grida)):
                gross_savings[age,z,a] = grida[int(choice_a[age,z,a])] - grida[a]
                savings_rate[age,z,a] = gross_savings[age,z,a] / (gridz[z,age]+r*grida[a])
    
    return gross_savings, savings_rate

# REWRITE THIS TO USE ARBITRARY QUANTS
def plot_savings_rate(age_start, n, gross_savings, savings_rate, distr_mass, mass_z):
    
    avg_rate_by_age_z = np.sum(savings_rate * distr_mass[:-1,:,:], axis = 2) / np.sum(distr_mass[:-1,:,:], axis = 2)
    avg_gross_by_age_z = np.sum(gross_savings * distr_mass[:-1,:,:], axis = 2) / np.sum(distr_mass[:-1,:,:], axis = 2)
    
    z_cum = np.cumsum(mass_z)
    
    # Plotting
    age_tick = np.arange(age_start,age_start+n)
    
    fig, (ax1,ax2) = plt.subplots(1,2, sharex = True, figsize=(12,6))
    for z in range(len(z_cum)):
        ax1.plot(age_tick, avg_rate_by_age_z[:,z], label = "Income quantile "+str(z_cum[z]))
        ax2.plot(age_tick, avg_gross_by_age_z[:,z], label = "Income quantile "+str(z_cum[z]))
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper left')
    ax2.set_xlabel('Household head age')
    ax1.set_ylabel('Average savings rate')
    ax2.set_ylabel('Average gross savings [BRL]')
    ax2.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.show()
    