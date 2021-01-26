
# Temptation_Inequality

<!-- badges: start -->
<!-- badges: end -->

Author: Felipe Bordini do Amaral

This repository holds the files used for my Masters' thesis

1. R/POF/tidying_data.R file Tidying POF data

  Requires .RDS files from "Pesquisa de Orçamento Familiares 2018"
  Requires 'Cadastro de Produto.xls' file, which categorizes all consumption as temptation or non-temptation
  
  Returns 'consumption_total.csv' containing item-by-item annualized consumption by household id
  Returns 'morador_count.csv' containing number of members by household id

2. .R file Tidying PNAD data and creating transitions and state-spaces

  Requires 'consumption_total.csv' and 'morador_count.csv' files from step 1.
  
  Returns 'temptation_filtered.csv' containing total, temptation and non-temptation consumption per capita by household id, after making pertinent data-cleaning
  
  Runs regression for case A and C, requires input from GMM for case B
  Plots comparison graphs for the three calibrations

3. .R file conducting parametric estimates


4. .py file with heterogeneous-agent model, calibration and graphs
