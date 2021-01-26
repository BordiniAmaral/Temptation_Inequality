# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 12:40:59 2020

@author: felip
"""
#---------------------------------- DEBUGGING CALIBRATION -----------------------------------------

import numpy as np
import CalibrationB as cb

pof_path = "C:\Backup Acer\EESP-Mestrado\Dissertação\POF\R Project\Exports\\temptation_filtered.csv"
pof_df = cb.import_POF_data(pof_path)
data_x, data_y = cb.select_data(pof_df)

x0_1 = 380*12 
sigx_1 = 3.09
sigy_1 = 1.92
xi_1 = 2.89

x0_2 = 300*12 
sigx_2 = 3.48
sigy_2 = 1.7
xi_2 = 2.7

gridq = np.arange(0,1,0.1)
W = np.identity(len(gridq))

# Running for calibration 1

consumpt, cons_bin_1, y_data_sel_sum_1, gridm_1 = cb.initial_computations(data_x, data_y, gridq, np.array([x0_1]))
y_sim_1 = np.zeros(shape = (len(data_x)))
x_sim_1 = np.zeros(shape = (len(data_x)))

for x in range(len(data_x)):
    if consumpt[x] < (x0_1 + 2):
        x_sim_1[x] = consumpt[x]
    else:
        x_sim_1[x] = cb.Newton(cb.f, cb.df, (x0_1 + 0.0001), 1e-6,int(1e6), consumpt[x], xi_1, sigx_1, sigy_1, x0_1)
        y_sim_1[x] = (data_x[x] + data_y[x]) - x_sim_1[x]

m_errors_1 = np.zeros(shape = len(gridm_1))
y_sim_sel_sum_1 = np.zeros(shape = len(gridm_1))
cons_bin_selected_1 = cons_bin_1[:,0]

for b in range(len(m_errors_1)):
    y_sim_sel_sum_1[b] = np.sum(y_sim_1[cons_bin_selected_1 == b])

m_errors_1 = ((y_sim_sel_sum_1 - y_data_sel_sum_1[:,0]) / y_data_sel_sum_1[:,0])
sqe_1 = (m_errors_1.reshape((1,len(gridm_1))) @ W @ m_errors_1.reshape((len(gridm_1),1)))[0][0]

# Running for calibration 2

consumpt, cons_bin_2, y_data_sel_sum_2, gridm_2 = cb.initial_computations(data_x, data_y, gridq, np.array([x0_2]))
y_sim_2 = np.zeros(shape = (len(data_x)))
x_sim_2 = np.zeros(shape = (len(data_x)))

for x in range(len(data_x)):
    if consumpt[x] < (x0_2 + 2):
        x_sim_2[x] = consumpt[x]
    else:
        x_sim_2[x] = cb.Newton(cb.f, cb.df, (x0_2 + 0.0001), 1e-6,int(1e6), consumpt[x], xi_2, sigx_2, sigy_2, x0_2)
        y_sim_2[x] = (data_x[x] + data_y[x]) - x_sim_2[x]

m_errors_2 = np.zeros(shape = len(gridm_2))
y_sim_sel_sum_2 = np.zeros(shape = len(gridm_2))
cons_bin_selected_2 = cons_bin_2[:,0]

for b in range(len(m_errors_2)):
    y_sim_sel_sum_2[b] = np.sum(y_sim_2[cons_bin_selected_2 == b])

m_errors_2 = ((y_sim_sel_sum_2 - y_data_sel_sum_2[:,0]) / y_data_sel_sum_2[:,0])
sqe_2 = (m_errors_2.reshape((1,len(gridm_2))) @ W @ m_errors_2.reshape((len(gridm_2),1)))[0][0]

# Comparing both at same time

sqe_both = cb.run_x0_simulations(np.array([sigx_1,sigx_2]), np.array([sigy_1,sigy_2]), np.array([x0_1,x0_2]), data_x, data_y, gridq, np.array([xi_1,xi_2]), W)
sqe_both_2 = np.array([sqe_both[0,0,0,0],sqe_both[1,1,1,1]])

sol = np.argwhere(sqe_both == np.min(sqe_both))
sol_2 = np.argwhere(sqe_both_2 == np.min(sqe_both_2))

sol_x0_1 = np.array([x0_1,x0_2])[sol[0][0]]
sol_sigx_1 = np.array([sigx_1,sigx_2])[sol[0][1]]
sol_sigy_1 = np.array([sigy_1,sigy_2])[sol[0][2]]
sol_xi_1 = np.array([xi_1,xi_2])[sol[0][3]]

# Oh well... running first stage entirely to see WHERE it fucking goes wrong - I guess I fixed it

grid_sigx = np.arange(3.45,3.9,0.01)
grid_sigy = np.arange(1.8,2.2,0.01)
grid_x0 = np.arange(250*12,310*12,5*12)
grid_xi = np.arange(3.1,3.5,0.01)

W = np.identity(len(gridq))

sqe_I = cb.run_x0_simulations(grid_sigx, grid_sigy, grid_x0, data_x, data_y, gridq, grid_xi, W)

sol = np.argwhere(sqe_I == np.min(sqe_I))
sol_x0_1 = grid_x0[sol[0][0]]
sol_sigx_1 = grid_sigx[sol[0][1]]
sol_sigy_1 = grid_sigy[sol[0][2]]
sol_xi_1 = grid_xi[sol[0][3]]

#---------------------------------- DEBUGGING OLG -----------------------------------------


