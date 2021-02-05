#-----------------------------------------------------------------------
# This file is dedicated to reading and tidying from Pesquisa de
# Orçamentos Familiares (POF) 2017-2018
#
# [REQUIRES]
# POF microdata compiled into RDS files following the script "Leitura dos Microdados - R.R"
#
# File by: Felipe Bordini do Amaral
# Institution: FGV - EESP
# Project: Master Dissertation
#
#-----------------------------------------------------------------------
#                     STARTING NECESSARY PACKAGES
#-----------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(glue)
library(nleqslv)
library(readxl)
library(stargazer)
library(cowplot)

#-----------------------------------------------------------------------
#                     DEFINING DATA PATH
#-----------------------------------------------------------------------

# Defining POF data path
dados_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/POF/2018/Dados_20200403"

# Defining export path for generated files
export_path <- "D:/Google Drive/Working Cloud//EESP-Mestrado/Dissertação/POF/R Project/Exports"

#-----------------------------------------------------------------------
#                     READING RDS FILES
#-----------------------------------------------------------------------

# Fetching rental tables
aluguel_rds <- readRDS(glue("{dados_path}/ALUGUEL_ESTIMADO.rds"))

# Fetching caderneta coletiva
caderneta_rds <- readRDS(glue("{dados_path}/CADERNETA_COLETIVA.rds"))

# Fetching collective expenditure
despesa_col_rds <- readRDS(glue("{dados_path}/DESPESA_COLETIVA.rds"))

# Fetching individual expenditure
despesa_ind_rds <- readRDS(glue("{dados_path}/DESPESA_INDIVIDUAL.rds"))

# Fetching household data
domicilio_rds <- readRDS(glue("{dados_path}/DOMICILIO.rds"))

# Fetching tenant data
morador_rds <- readRDS(glue("{dados_path}/MORADOR.rds"))

# Fetching other incomes
outros_rendimentos_rds <- readRDS(glue("{dados_path}/OUTROS_RENDIMENTOS.rds"))

# Fetching main incomes
rendimentos_rds <- readRDS(glue("{dados_path}/RENDIMENTO_TRABALHO.rds"))

# Fetching consumption/income variable library
# This file contains the original labeling of products/incomes, plus two columns:
# - User-defined categories
# - User-defined temptation status

cadastro_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/POF/2018/Documentacao_20200403/Cadastro de Produtos.xls"
cadastro_produto <- read_xls(cadastro_path, col_types = "text")

temp_cadastro <- cadastro_produto %>% select(V9001,V9000,category,temptation)

#-----------------------------------------------------------------------
#                     ADDING CATEGORY/TEMPTATION TO DATABASE
#-----------------------------------------------------------------------

# Adding category/temptation to expenditure tables
despesa_ind_aggr <- despesa_ind_rds %>%
  mutate(V9001 = as.character(V9001)) %>%
  left_join(temp_cadastro, by = "V9001")

despesa_col_aggr <- despesa_col_rds %>%
  mutate(V9001 = as.character(V9001)) %>%
  left_join(temp_cadastro, by = "V9001")

caderneta_aggr <- caderneta_rds %>%
  mutate(V9001 = as.character(V9001)) %>%
  left_join(temp_cadastro, by = "V9001")

aluguel_aggr <- aluguel_rds %>%
  mutate(V9001 = str_pad(V9001, 6, pad = "0")) %>%
  left_join(temp_cadastro, by = "V9001")

selling_aggr <- outros_rendimentos_rds %>%
  mutate(V9001 = as.character(V9001)) %>%
  left_join(temp_cadastro, by = "V9001") %>%
  filter(category == "durable selling") %>%
  mutate(V8501_DEFLA = if_else(is.na(V8501_DEFLA),0,V8501_DEFLA),
         transaction = -(V8500_DEFLA - V8501_DEFLA),
         temptation = "non-temptation")


#-----------------------------------------------------------------------
#                     RETRIEVING HOUSEHOLD INFO
#-----------------------------------------------------------------------

# Counting how many people per household

morador_aggr <- morador_rds %>%
  mutate(COD_UPA = str_pad(COD_UPA,9,pad = "0"),
         NUM_DOM = str_pad(NUM_DOM,2,pad = "0"),
         NUM_UC = str_pad(NUM_UC,2,pad = "0"),
         hh_id = paste0(COD_UPA,NUM_DOM,NUM_UC))

morador_count <- morador_aggr %>%
  group_by(hh_id) %>%
  summarise(n = n())

#-----------------------------------------------------------------------
#                     JOINING DATABASES
#-----------------------------------------------------------------------

cad_select <- caderneta_aggr %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,V8000_DEFLA,FATOR_ANUALIZACAO,V9000,category,temptation) %>%
  rename(transaction = V8000_DEFLA)
col_select <- despesa_col_aggr %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,V8000_DEFLA,FATOR_ANUALIZACAO,V9000,category,temptation) %>%
  rename(transaction = V8000_DEFLA)
ind_select <- despesa_ind_aggr %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,V8000_DEFLA,FATOR_ANUALIZACAO,V9000,category,temptation) %>%
  rename(transaction = V8000_DEFLA)
sel_select <- selling_aggr %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,transaction,FATOR_ANUALIZACAO,V9000,category,temptation)
alg_select <- aluguel_aggr %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,V8000_DEFLA,FATOR_ANUALIZACAO,V9000,category,temptation) %>%
  rename(transaction = V8000_DEFLA)

# Setting up appropriate names to keep track where each row comes from
consumption_total <- bind_rows(cad_select,
                               col_select,
                               ind_select,
                               sel_select,
                               alg_select,
                               .id = "origin")
names <- c("caderneta","coletiva","individual","selling","aluguel")

consumption_total <- consumption_total %>%
  mutate(origin = names[as.integer(origin)])

#-----------------------------------------------------------------------
#                     AGGREGATING BY ID
#-----------------------------------------------------------------------

# aggregating all consumption into temptation/non-temptation by Consumption Unit
# [1]: COD_UPA
# [2]: NUM_DOM
# [3]: NUM_UC
# Unique identifier: 11112233

consumption_total <- consumption_total %>%
  mutate(COD_UPA = str_pad(COD_UPA,9,pad = "0"),
         NUM_DOM = str_pad(NUM_DOM,2,pad = "0"),
         NUM_UC = str_pad(NUM_UC,2,pad = "0"),
         id = paste0(COD_UPA,NUM_DOM,NUM_UC),
         annualized_exp = transaction * FATOR_ANUALIZACAO)

# consumption_total is the main source of data for all exercises

#-----------------------------------------------------------------------
#                     INCOME TABLES
#-----------------------------------------------------------------------

cadastro_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/POF/2018/Documentacao_20200403/Cadastro de Produtos.xls"
cadastro_produto <- read_xls(cadastro_path, col_types = "text")
temp_cadastro <- cadastro_produto %>% select(V9001,V9000,category)

# Adding category to "other incomes"
outros_rendimentos_rds <- outros_rendimentos_rds %>%
  mutate(V9001 = as.character(V9001)) %>%
  left_join(temp_cadastro, by = "V9001")

# Selecting pertinent columns from income tables

rendimento_select <- rendimentos_rds %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,V9011,V8500_DEFLA,V531112_DEFLA,V531122_DEFLA,V531132_DEFLA) %>%
  replace_na(list(V9011 = 1)) %>%
  replace(is.na(.),0) %>%
  mutate(category = "income",
         COD_UPA = str_pad(COD_UPA,9,pad = "0"),
         NUM_DOM = str_pad(NUM_DOM,2,pad = "0"),
         NUM_UC = str_pad(NUM_UC,2,pad = "0"),
         id = paste0(COD_UPA,NUM_DOM,NUM_UC),
         annual_defla = V9011*(V8500_DEFLA - V531112_DEFLA - V531122_DEFLA - V531132_DEFLA))

outros_rendimentos_select <- outros_rendimentos_rds %>%
  select(COD_UPA,NUM_DOM,NUM_UC,V9001,V9011,V8500_DEFLA,V8501_DEFLA, category) %>%
  replace_na(list(V9011 = 1)) %>%
  replace(is.na(.),0) %>%
  mutate(COD_UPA = str_pad(COD_UPA,9,pad = "0"),
         NUM_DOM = str_pad(NUM_DOM,2,pad = "0"),
         NUM_UC = str_pad(NUM_UC,2,pad = "0"),
         id = paste0(COD_UPA,NUM_DOM,NUM_UC),
         annual_defla = V9011*(V8500_DEFLA - V8501_DEFLA))

# Grouping by household all relevant incomes and expenditures in order to obtain
# savings by Income - Consumption

outros_rendimentos_clean <- outros_rendimentos_select %>%
  filter(category %in% c("income","bequest","transfers")) %>%
  select(id,annual_defla)

rendimento_clean <- rendimento_select %>%
  select(id,annual_defla)

all_income <- rendimento_clean %>%
  bind_rows(outros_rendimentos_clean) %>%
  group_by(id) %>%
  summarise(total_income = sum(annual_defla)) %>%
  ungroup() %>%
  left_join(morador_count, by = c("id" = "hh_id")) %>%
  mutate(income_pc = total_income / n)

all_consumption <- consumption_total %>%
  filter(temptation != "n/a") %>%
  filter(!(category %in% c("household acquisition",
                           "durable selling",
                           "vehicle acquisition",
                           "household construction and repair",
                           "social security"))) %>%
  select(id,annualized_exp) %>%
  group_by(id) %>%
  summarise(total_expenditure = sum(annualized_exp)) %>%
  ungroup() %>%
  left_join(morador_count, by = c("id" = "hh_id")) %>%
  mutate(exp_pc = total_expenditure / n)

savings_by_residual <- all_consumption %>%
  select(-n) %>%
  inner_join(all_income, by = "id") %>%
  mutate(gross_savings_pc = income_pc - exp_pc,
         savings_rate = gross_savings_pc / income_pc)

#-----------------------------------------------------------------------
#                     SAVING TABLES
#-----------------------------------------------------------------------

 setwd(export_path)
 fwrite(consumption_total, "consumption_total.csv", sep = ";")
 fwrite(savings_by_residual, "savings_by_residual.csv", sep = ";")
 fwrite(morador_count, "morador_count.csv", sep = ";")
