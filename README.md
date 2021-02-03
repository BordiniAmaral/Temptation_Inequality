
# Temptation_Inequality

<!-- badges: start -->
<!-- badges: end -->

Author: Felipe Bordini do Amaral

This repository holds the files used for my Masters' thesis

## 1. R/POF/tidying_data.R file Tidying POF data

  Requires .RDS files from "Pesquisa de Orçamento Familiares 2018"
  Requires 'Cadastro de Produto.xls' file, which categorizes all consumption as temptation or non-temptation
  
  Returns 'consumption_total.csv' containing item-by-item annualized consumption by household id
  Returns 'morador_count.csv' containing number of members by household id

## 2. R/POF/Stoney-Geary.R file organizing and comparing calibrations A, B and C

  Requires 'consumption_total.csv' and 'morador_count.csv' files from step 1.
  
  Returns 'temptation_filtered.csv' containing total, temptation and non-temptation consumption per capita by household id, after making pertinent data-cleaning
  
  Runs regression for case A and C, requires input from GMM for case B
  Plots comparison graphs for the three calibrations

## 3. R/POF/Graphics.R file generating additional plots

  Requires 'consumption_total.csv' and 'morador_count.csv' files from step 1.


## 4. R/PNAD/Main.R file Tidying PNAD data and creating transitions and state-spaces

  Imports and saves PNAD files (do only once, it takes some time)
  Requires 'deflator_PNADC_2020_trimestral_010203.xls' file from IBGE
  Requires 'TributosPrevidencia.xlsx' with tax and social security schedules
  Deducts individual tax and social security schedule
  
  Returns 'quantile_frontiers.csv' and 'quantile_values.csv' defining the measures and values for the disposable income per capita state-space
  Returns 'complete_transitions.csv' with all PNADC (t,t+4) transitions for the period
  Returns matrices and graphs for transition (by age and aggregate)

## 5. Python/ .py files with heterogeneous-agent model, calibration and graphs:

  Main.py is the main file, calling all others
  Calibration.py conducts a blunt force GMM calibration of parameters x0, gamma1 and gamma2
  Graphics.py executes some auxiliary calculations and plot resulting equilibrium stats
  OLG_A_CRRA.py solves the lifecycle model in general equilibrium
  Temptation_CRRA.py solves allocations for temptation and non-temptation consumptions
  Income_Process.py imports the .csv file with transition matrix and income values
