#-----------------------------------------------------------------------
# This file is dedicated to building the income transition matrices for
# Brasil, based upon PNADc data.
#
#
# File by: Felipe Bordini do Amaral
# Institution: FGV - EESP
# Project: Master Dissertation
#
#-----------------------------------------------------------------------
#                     STARTING NECESSARY PACKAGES
#-----------------------------------------------------------------------

library(tidyverse)
library(PNADcIBGE)
library(readxl)
library(glue)
library(data.table)

#-----------------------------------------------------------------------
#                     DOWNLOADING AND TIDYING PNADc DATA
#-----------------------------------------------------------------------


# Defining tax/social security schedules discounts

IR_discount <- function(public, formal, value, schedule){
  return(case_when(public == 1 | formal == 1 ~ case_when(value > max(schedule$ref_upper) ~ value * schedule$rate[which.max(schedule$lower)] - schedule$deduction[which.max(schedule$lower)],
                                                            value < min(schedule$lower) ~ 0,
                                                            TRUE ~ ifelse(is_empty(which(schedule$lower < value & schedule$upper > value)), -1, value * schedule$rate[which(schedule$lower < value & schedule$upper > value)] - schedule$deduction[which(schedule$lower < value & schedule$upper > value)])),
                      TRUE ~ 0))}

INSS_discount <- function(public, formal, contributes, value, schedule){
  return(case_when(public == 1 ~ value * 0.2,
                   formal == 1 | contributes == 1 ~ case_when(value > max(schedule$ref_upper) ~ max(schedule$ref_upper) * schedule$rate[which.max(schedule$lower)],
                                                              TRUE ~ ifelse(is_empty(which(schedule$lower < value & schedule$upper > value)), -1, value * schedule$rate[which(schedule$lower < value & schedule$upper > value)])),
                   TRUE ~ 0))}

# Adpating first observation for the loop
# Create a 13-digit id for household joining:
## UPA (9 digits) + V1008 (2 digits) + V1014 (2 digits)
# Selecting number of people by household (V2001)
# Selecting age of household head (V2009 for V2005 == 01)
# Adding respective tax schedule
# aggregate household income net of taxes
## V403412 - actual gross income in main job
## V403422 - actual gross good income in main job
## V405112 - actual gross income in secondary job
## V405122 - actual gross good income in secondary job
## V405912 - actual gross income in other jobs
## V405922 - actual gross good income in other jobs
### Variables for determining disposable income:
### V4028 - 1 if main job was in public sector
### V4029 - 1 if main job was formal
### V4032 - 1 if main job contributed to social security
### V4047 - 1 if secondary job was in public sector
### V4048 - 1 if secondary job was formal
### V4049 - 1 if secondary job contributed to social security
### V4057 - 1 if other jobs contributed to social security

# Defining function for processing PNADc quarterly data

pnadc_income_periods <- function(periods,tax_table,export_path){

  # Checking order of periods

  periods <- periods %>% arrange(ano,tri)

  for(t in 1:length(periods$tri)){

    cat(glue("------------------------------ \n Starting year {periods$ano[t]}, quarter {periods$tri[t]} \n \n"))

    # Selecting tax schedule for given period
    IR_temp <- tax_table %>%
      filter(ano == periods$ano[t],
             tri == periods$tri[t] | tri == "all",
             type == "IR")

    INSS_temp <- tax_table %>%
      filter(ano == periods$ano[t],
             tri == periods$tri[t] | tri == "all",
             type == "INSS")

    cat("Fetching PNADc observation... \n")
    pnadc_df <- get_pnadc(year = periods$ano[t],
                          quarter = periods$tri[t],
                          design = FALSE,
                          labels = FALSE)

    current_defla <- pnadc_df %>%
      select(UF, Efetivo) %>%
      distinct()

    price_deflator <- ref_prices %>%
      left_join(current_defla, by = "UF") %>%
      mutate(defla = Efetivo / EfetivoRef) %>%
      select(UF, defla)

    cat("Working on PNADC observation... \n")
    current_df <- pnadc_df %>%
      mutate(hh_id = str_c(str_pad(UPA,9,pad = "0"), str_pad(V1008,2,pad = "0"), str_pad(V1014,9,pad = "0"))) %>%
      select(Ano, Trimestre, UF, hh_id, V1016,
             VD2003, V2005, V2009,
             V4028, V4029, V4032, V403412, V403422,
             V4047, V4048, V4049, V405112, V405122,
             V4057, V405912, V405922) %>%
      left_join(price_deflator, by = "UF") %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      rowwise() %>%
      mutate(main_IR = IR_discount(V4028,V4029,V403412,IR_temp),
             sec_IR =  IR_discount(V4047,V4048,V405112,IR_temp),
             main_INSS = INSS_discount(V4028,V4029,V4032,V403412,INSS_temp),
             sec_INSS = INSS_discount(V4047,V4048,V4049,V405112,INSS_temp),
             other_INSS = INSS_discount(0,0,V4057,V405912,INSS_temp),
             total_income = (V403412+V403422+V405112+V405122+V405912+V405922)*12 + 1,33*(V403412*(V4029 == 1) + V405112*(V4048 == 1)),
             total_deduction = (main_IR + sec_IR + main_INSS + sec_INSS + other_INSS)*12,
             total_net = total_income - total_deduction,
             total_net_defla = total_net * defla,
             head = (V2005 == "01")) %>%
      ungroup() %>%
      group_by(hh_id) %>%
      summarise(hh_net_defla_annualized = sum(total_net_defla),
                hh_age_head = max(head*V2009),
                hh_members = max(VD2003),
                V1016 = max(V1016)) %>%
      mutate(hh_net_defla_annualized_pc = hh_net_defla_annualized / hh_members,
             t = t)

    cat("Saving current observation... \n")

    setwd(export_path)
    fwrite(current_df, file = glue("PNADC_clean_{periods$ano[t]}_{periods$tri[t]}.csv"), sep = ";")
  }
}

# Using jan/2018 prices to match POF 2017-2018 prices
# Getting deflators for 1Q2018:
price_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/deflator_PNADC_2020_trimestral_010203.xls"
ref_prices <- read_xls(price_path) %>%
  filter(Ano == 2018,
         trim == "01-02-03") %>%
  select(UF, Efetivo) %>%
  rename(EfetivoRef = Efetivo)
rm(price_path)

# Getting tax and social security contribution schedules:
tax_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/TributosPrevidencia.xlsx"
tax_table <- read_xlsx(tax_path) %>%
  replace_na(list(upper = Inf, deduction = 0))
rm(tax_path)

periods <- data.frame(ano = rep(2012:2018, each = 4),
                      tri = rep(1:4, (2018-2012)+1))

# Running PNADc observations
export_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/R Project/Clean household income"
pnadc_income_periods(periods,tax_table,export_path)

#-----------------------------------------------------------------------
#                     JOINING ALL OBSERVATIONS
#-----------------------------------------------------------------------

# Imports all files from selected folder
# Make sure only PNADC_clean csv files are there

import_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/R Project/Clean household income"
import_files <- list.files(path = import_path)

pnadc_hh <- map(glue("{import_path}/{import_files}"), fread, sep = ";", dec = ".")
pnadc_hh <- bind_rows(pnadc_hh)

#-----------------------------------------------------------------------
#      Establishing state-space for quantiles of hh income by age
#-----------------------------------------------------------------------

export_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/R Project/Transitions by age"
setwd(export_path)

quantile_frontiers <- pnadc_hh %>%
  filter(hh_age_head >= 25, hh_age_head <= 66, V1016 %in% c(1,5)) %>%
  group_by(hh_id) %>%
  mutate(both = (n()==2),
         same_ref = (max(hh_age_head)-min(hh_age_head)) <= 1) %>%
  ungroup() %>%
  filter(both == TRUE,
         same_ref == TRUE) %>%
  group_by(hh_age_head) %>%
  summarise(q98 = quantile(hh_net_defla_annualized_pc, probs = 0.98, type = 8),
            q96 = quantile(hh_net_defla_annualized_pc, probs = 0.96, type = 8),
            q94 = quantile(hh_net_defla_annualized_pc, probs = 0.94, type = 8),
            q92 = quantile(hh_net_defla_annualized_pc, probs = 0.92, type = 8),
            q90 = quantile(hh_net_defla_annualized_pc, probs = 0.90, type = 8),
            q80 = quantile(hh_net_defla_annualized_pc, probs = 0.80, type = 8),
            q70 = quantile(hh_net_defla_annualized_pc, probs = 0.70, type = 8),
            q60 = quantile(hh_net_defla_annualized_pc, probs = 0.60, type = 8),
            q50 = quantile(hh_net_defla_annualized_pc, probs = 0.50, type = 8),
            q40 = quantile(hh_net_defla_annualized_pc, probs = 0.40, type = 8),
            q30 = quantile(hh_net_defla_annualized_pc, probs = 0.30, type = 8),
            q20 = quantile(hh_net_defla_annualized_pc, probs = 0.20, type = 8),
            q10 = quantile(hh_net_defla_annualized_pc, probs = 0.10, type = 8)) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quant")

quantile_values <- pnadc_hh %>%
  filter(hh_age_head >= 25, hh_age_head <= 66, V1016 %in% c(1,5)) %>%
  group_by(hh_id) %>%
  mutate(both = (n()==2),
         same_ref = (max(hh_age_head)-min(hh_age_head)) <= 1) %>%
  ungroup() %>%
  filter(both == TRUE,
         same_ref == TRUE) %>%
  group_by(hh_age_head) %>%
  summarise(v98 = quantile(hh_net_defla_annualized_pc, probs = 0.99, type = 8),
            v96 = quantile(hh_net_defla_annualized_pc, probs = 0.97, type = 8),
            v94 = quantile(hh_net_defla_annualized_pc, probs = 0.95, type = 8),
            v92 = quantile(hh_net_defla_annualized_pc, probs = 0.93, type = 8),
            v90 = quantile(hh_net_defla_annualized_pc, probs = 0.91, type = 8),
            v80 = quantile(hh_net_defla_annualized_pc, probs = 0.85, type = 8),
            v70 = quantile(hh_net_defla_annualized_pc, probs = 0.75, type = 8),
            v60 = quantile(hh_net_defla_annualized_pc, probs = 0.65, type = 8),
            v50 = quantile(hh_net_defla_annualized_pc, probs = 0.55, type = 8),
            v40 = quantile(hh_net_defla_annualized_pc, probs = 0.45, type = 8),
            v30 = quantile(hh_net_defla_annualized_pc, probs = 0.35, type = 8),
            v20 = quantile(hh_net_defla_annualized_pc, probs = 0.25, type = 8),
            v10 = quantile(hh_net_defla_annualized_pc, probs = 0.15, type = 8),
            v00 = quantile(hh_net_defla_annualized_pc, probs = 0.05, type = 8)) %>%
  pivot_longer(cols = starts_with("v"), names_to = "quant_v")

fwrite(quantile_frontiers, file = "quantile_frontiers.csv", sep = ";")
fwrite(quantile_values, file = "quantile_values.csv", sep = ";")

# Quantile frontiers current (t) and next (t+4) ####

quantile_frontiers_current <- pnadc_hh %>%
  filter(V1016 %in% c(1,5)) %>%
  group_by(hh_id) %>%
  mutate(both = (n()==2)) %>%
  ungroup() %>%
  filter(both == TRUE,
         V1016 == 1) %>%
  filter(hh_age_head >= 25, hh_age_head <= 65) %>%
  group_by(hh_age_head) %>%
  summarise(q98 = quantile(hh_net_defla_annualized_pc, probs = 0.98, type = 8),
            q96 = quantile(hh_net_defla_annualized_pc, probs = 0.96, type = 8),
            q94 = quantile(hh_net_defla_annualized_pc, probs = 0.94, type = 8),
            q92 = quantile(hh_net_defla_annualized_pc, probs = 0.92, type = 8),
            q90 = quantile(hh_net_defla_annualized_pc, probs = 0.90, type = 8),
            q80 = quantile(hh_net_defla_annualized_pc, probs = 0.80, type = 8),
            q70 = quantile(hh_net_defla_annualized_pc, probs = 0.70, type = 8),
            q60 = quantile(hh_net_defla_annualized_pc, probs = 0.60, type = 8),
            q50 = quantile(hh_net_defla_annualized_pc, probs = 0.50, type = 8),
            q40 = quantile(hh_net_defla_annualized_pc, probs = 0.40, type = 8),
            q30 = quantile(hh_net_defla_annualized_pc, probs = 0.30, type = 8),
            q20 = quantile(hh_net_defla_annualized_pc, probs = 0.20, type = 8),
            q10 = quantile(hh_net_defla_annualized_pc, probs = 0.10, type = 8)) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quant")

quantile_frontiers_next <- pnadc_hh %>%
  filter(V1016 %in% c(1,5)) %>%
  group_by(hh_id) %>%
  mutate(both = (n()==2),
         age_1 = max((V1016 == 1)*hh_age_head)) %>%
  ungroup() %>%
  filter(both == TRUE,
         V1016 == 5) %>%
  mutate(hh_age_head = age_1) %>%
  select(-age_1) %>%
  filter(hh_age_head >= 25, hh_age_head <= 65) %>%
  group_by(hh_age_head) %>%
  summarise(q98 = quantile(hh_net_defla_annualized_pc, probs = 0.98, type = 8),
            q96 = quantile(hh_net_defla_annualized_pc, probs = 0.96, type = 8),
            q94 = quantile(hh_net_defla_annualized_pc, probs = 0.94, type = 8),
            q92 = quantile(hh_net_defla_annualized_pc, probs = 0.92, type = 8),
            q90 = quantile(hh_net_defla_annualized_pc, probs = 0.90, type = 8),
            q80 = quantile(hh_net_defla_annualized_pc, probs = 0.80, type = 8),
            q70 = quantile(hh_net_defla_annualized_pc, probs = 0.70, type = 8),
            q60 = quantile(hh_net_defla_annualized_pc, probs = 0.60, type = 8),
            q50 = quantile(hh_net_defla_annualized_pc, probs = 0.50, type = 8),
            q40 = quantile(hh_net_defla_annualized_pc, probs = 0.40, type = 8),
            q30 = quantile(hh_net_defla_annualized_pc, probs = 0.30, type = 8),
            q20 = quantile(hh_net_defla_annualized_pc, probs = 0.20, type = 8),
            q10 = quantile(hh_net_defla_annualized_pc, probs = 0.10, type = 8)) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quant")

#-----------------------------------------------------------------------
#                    Creating (t, t+4) pairing
#-----------------------------------------------------------------------

export_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/R Project/Transitions by age"
setwd(export_path)

compile_pairwise_transition <- function(pnadc_hh, periods){

  # Starting final dataframe

  complete_df <- data.frame(hh_id = character(),
                            hh_members = integer(),
                            hh_age_head = integer(),
                            y_current = double(),
                            y_next = double())

  cat("-------------------------------------\n")
  for(i in 1:periods){

    if(i+4>periods){break}

    cat(glue("  Running period {i} out of {periods-4} \n \n"))
    # Does not run at pairwise start-up

    current_df <- pnadc_hh %>%
      filter(t == i, V1016 == 1) %>%
      rename(y_current = hh_net_defla_annualized_pc) %>%
      select(hh_id, hh_age_head, hh_members, y_current)

    next_df <- pnadc_hh %>%
      filter(t == i+4, V1016 == 5) %>%
      rename(y_next = hh_net_defla_annualized_pc) %>%
      select(hh_id, y_next)

    joined_df <- current_df %>%
      left_join(next_df, by = "hh_id") %>%
      filter(!is.na(y_next))

    # Adds to complete data-frame

    complete_df <- bind_rows(complete_df, joined_df)

  }
  cat("             Finish \n")
  cat("-------------------------------------")
  return(complete_df)
}

complete_transitions <- compile_pairwise_transition(pnadc_hh, periods = 28)
fwrite(complete_transitions, file = "complete_transitions.csv")

#-----------------------------------------------------------------------
#            Calculating implied transition matrices by age
#-----------------------------------------------------------------------

ages <- seq(from = 25, to = 65, by = 1)
export_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/PNADC/R Project/Transitions by age"

build_transitions <- function(age,
                              quantile_frontiers,
                              quantile_values,
                              complete_transitions,
                              export_path){

  cat("----------------------\n")
  cat(glue("Starting for age {age} \n"))

  setwd(export_path)

  current_frontier <- quantile_frontiers_current %>%
    filter(hh_age_head == age) %>%
    arrange(value)

  next_frontier <- quantile_frontiers_next %>%
    filter(hh_age_head == age) %>%
    arrange(value)

  current_transition <- complete_transitions %>%
    select(-c("hh_id")) %>%
    filter(hh_age_head == age) %>%
    rowwise() %>%
    mutate(current_q = sum(current_frontier$value < y_current),
           next_q = sum(next_frontier$value < y_next))

  transition_matrix <- current_transition %>%
    count(current_q, next_q) %>%
    group_by(current_q) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = next_q, names_prefix = "q", values_from = freq, values_fill = 0) %>%
    arrange(current_q)

  # Building a proper matrix for this transition
  row.names(transition_matrix) <- glue("q{transition_matrix$current_q}")
  quantiles <- seq(from = 0, to = 13, by = 1)
  matrix_1 <- as.matrix(transition_matrix[,glue("q{quantiles}")])

  fwrite(matrix_1, file = glue("{age}yo_matrix1.csv"))

  # Density plot of current transition

  quantile_points <- c(seq(0,0.9,0.1),
                       seq(0.92,1,0.02))

  current_transition <- current_transition %>%
    mutate(current_quantile = quantile_points[as.integer(current_q)+1],
           next_quantile = quantile_points[as.integer(next_q)+1])

  quantile_points <- data.frame(current_quantile = quantile_points,
                                x_vline = quantile_points)

  current_transition$current_quantile <- as.factor(current_transition$current_quantile)
  quantile_points$current_quantile <- as.factor(quantile_points$current_quantile)

  ggplot(current_transition)+
    geom_density(aes(next_quantile, colour = current_quantile, fill = current_quantile), adjust = 4, alpha = 0.1)+
    facet_wrap(~current_quantile) +
    geom_vline(data = quantile_points, aes(xintercept = x_vline), linetype = "dashed") +
    labs(title = glue("Quantile transition distributions - {age} y.o."),
         subtitle = "PNADC 2012-2018",
         x = "Next quantile",
         y = "Density",
         colour = "Current quantile",
         fill = "Current quantile") +
    ggsave(filename = glue("{age}yo_distribution.png"),
           width = 10,
           height = 8,
           dpi = 300)

  return(current_transition)
}

transitions_with_quantile <- map(ages,build_transitions,
                                 quantile_frontiers,
                                 quantile_values,
                                 complete_transitions,
                                 export_path) %>% bind_rows()

#-----------------------------------------------------------------------
#            Compiling an aggregate transition
#-----------------------------------------------------------------------
# joining all ages transitions ####

average_transition <- transitions_with_quantile %>%
  count(current_q, next_q) %>%
  group_by(current_q) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = next_q, names_prefix = "q", values_from = freq, values_fill = 0) %>%
  arrange(current_q)

# Building a proper matrix for the average transition ####
row.names(average_transition) <- glue("q{average_transition$current_q}")
quantiles <- seq(from = 0, to = 13, by = 1)
matrix_a <- as.matrix(average_transition[,glue("q{quantiles}")])

fwrite(matrix_a, file = glue("average_matrix1.csv"))

# Density plots of average transition ####

# Creating labels

quantile_points <- c(seq(0.1,0.9,0.1),
                     seq(0.92,1,0.02))

interval_l <- c(seq(0,0.9,0.1),
                seq(0.92,0.98,0.02))

interval_u <- c(seq(0.1,0.9,0.1),
                seq(0.92,1,0.02))

width <- interval_u - interval_l
quant.lab <- glue("{interval_l} - {interval_u}")
names(quant.lab) <- quantile_points

# Adding quantile intervals to dataframe

transitions_with_quantile <- transitions_with_quantile %>%
  mutate(current_quantile = quantile_points[as.integer(current_q)+1],
         next_quantile = quantile_points[as.integer(next_q)+1],
         mass = width[as.integer(next_q)+1])

quantile_points <- data.frame(current_quantile = quantile_points,
                              x_vline = quantile_points)

# Preparing for Plotting

transitions_with_quantile$current_quantile <- as.factor(transitions_with_quantile$current_quantile)
quantile_points$current_quantile <- as.factor(quantile_points$current_quantile)

# Plotting

ggplot(transitions_with_quantile)+
  geom_density(aes(next_quantile, colour = current_quantile, fill = current_quantile), kernel = "gaussian", adjust = 5, alpha = 0.1, outline.type = "both", trim = TRUE)+
  facet_wrap(~current_quantile, labeller = labeller(current_quantile = quant.lab), ncol = 5) +
  geom_vline(data = quantile_points, aes(xintercept = x_vline), linetype = "dashed") +
  labs(title = glue("Disposable income per capita one-year transition - Average 25-65 y.o. family heads"),
       subtitle = " PNADC 2012-2018 \n Facet label: Current Quantile",
       x = "Next quantile",
       y = "Density",
       colour = "Current quantile",
       fill = "Current quantile") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0,1,0.5)) +
  ggsave(filename = glue("IncomeDistr.png"),
         width = 8,
         height = 8,
         dpi = 600)

ggplot(transitions_with_quantile)+
  geom_density(aes(next_quantile, colour = current_quantile, fill = current_quantile), kernel = "gaussian", adjust = 5, alpha = 0.1, outline.type = "both", trim = TRUE)+
  facet_wrap(~current_quantile, labeller = labeller(current_quantile = quant.lab), ncol = 7) +
  geom_vline(data = quantile_points, aes(xintercept = x_vline), linetype = "dashed") +
  labs(title = glue("Disposable income per capita one-year transition - Average 25-65 y.o. family heads"),
       subtitle = " PNADC 2012-2018 \n Facet label: Current Quantile",
       x = "Next quantile",
       y = "Density",
       colour = "Current quantile",
       fill = "Current quantile") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0,1,0.5)) +
  ggsave(filename = glue("IncomeDistr7col.png"),
         width = 11,
         height = 6,
         dpi = 600)

# Compiling average # of household members by age and quantile

members_by_quantile <- transitions_with_quantile %>%
  group_by(hh_age_head, current_q) %>%
  summarise(mean_num = mean(hh_members))

setwd(export_path)
fwrite(members_by_quantile,"members_by_quantile_age.csv", sep = ";")

#-----------------------------------------------------------------------
#                     *** LE RABISQUERA ***
#-----------------------------------------------------------------------

# Checking the average transition matrix

quantile_mass <- c(rep(0.02,5),
                   rep(0.1,8),
                   rep(0.02,5))

deviation <- (quantile_mass %*% matrix_a) / quantile_mass

# Testing asymptotic distribution
a_new <- rep(1/18, 18)
error <- 1

while(error>1e-6){
  a <- a_new
  a_new <- a %*% matrix_a
  error <- max(abs(a_new-a))
}

# Checking transitions and #members by household

for_plot <- transitions_with_quantile
for_plot$current_q <- as.numeric(for_plot$current_q)

ggplot(for_plot)+
  geom_smooth(aes(x = hh_age_head, y = hh_members))+
  facet_wrap(~current_q)

for_plot <- members_by_quantile
for_plot$current_q <- as.factor(for_plot$current_q)

ggplot(for_plot)+
  geom_line(aes(x = hh_age_head, y = mean_num, colour = current_q))

# Checking distributions

pnadc_plot_selec <- pnadc_hh %>%
  filter(hh_age_head >= 25, hh_age_head <= 65)

pnadc_plot_selec$hh_age_head <- as.factor(pnadc_plot_selec$hh_age_head)

ggplot(pnadc_plot_selec) +
  geom_violin(aes(x = hh_age_head, y = hh_net_defla_annualized_pc/12), scale = "area")

