#-----------------------------------------------------------------------
# This file is dedicated to reading and tidying from Pesquisa de
# Orçamentos Familiares (POF) 2017-2018
#
# [REQUIRES] - 
# consumption_total dataframe obtained from tidying_data.R
# morador_count dataframe obtained from tidying_data.R
# 
# Recommended generate file from script instead of loading
# file to avoid export/import issues
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
#                     DEFINING EXPORTING PATH
#-----------------------------------------------------------------------

# Cleaning undesired objects
rm(list = setdiff(ls(), c("consumption_total","morador_count")))

export_path <- "C:/Backup Acer/EESP-Mestrado/Dissertação/POF/R Project/Exports"
setwd(export_path)

#-----------------------------------------------------------------------
#                     DEFINING EXERCISE FUNCTION
#-----------------------------------------------------------------------

regress_and_plot <- function(consumption_total,
                             morador_count,
                             ex_name,
                             trim_outliers = FALSE,
                             lower_trim = 0,
                             upper_trim = 100000,
                             remove_durables = FALSE,
                             save = FALSE){
  
  temptation_cross <- consumption_total %>%
    filter(temptation != "n/a")
  
  if(remove_durables){
    temptation_cross <-  temptation_cross %>%
      filter(!(category %in% c("household acquisition","durable selling","vehicle acquisition")))
  }
  
  temptation_cross <- temptation_cross %>%
    group_by(id, temptation) %>%
    summarise(expenditure = sum(annualized_exp)) %>%
    pivot_wider(id_cols = id, 
                names_from = temptation, 
                values_from = expenditure,
                values_fill = 0) %>%
    left_join(morador_count, by = c("id" = "hh_id")) %>%
    mutate(total = `non-temptation` + temptation,
           total_monthly = total/12,
           total_monthly_pc = total_monthly / n,
           temptation_pc = temptation /n,
           `non-temptation_pc` = `non-temptation` / n,
           tempt_frac = temptation/total)
  
  if(trim_outliers){
    temptation_cross <- temptation_cross %>%
      filter(total_monthly_pc > lower_trim, total_monthly_pc < upper_trim)
  }
  
  # plots
  
  plot1 <-  ggplot(temptation_cross) +
    geom_point(aes(x = total_monthly_pc, y = tempt_frac), size = 0.5) +
    labs(title = "Temptation by Consumption level",
         subtitle = "POF 2017 - 2018 (per-capita)",
         x = "Monthly expenditure per capita [R$]",
         y = "Fraction spent on temptation")
  
  
  plot2 <-  ggplot(temptation_cross) +
    geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
    labs(title = "Temptation by Consumption level",
         subtitle = "POF 2017 - 2018 (per-capita)",
         x = "Monthly expenditure per capita [R$]",
         y = "Fraction spent on temptation")
    
  
  plot3 <-  ggplot(temptation_cross) +
    geom_smooth(aes(x = log(`non-temptation_pc`), y = log(temptation_pc))) +
    labs(title = "Tempt. Cons. by Non-tempt. Cons.",
         subtitle = "POF 2017 - 2018 (per-capita)",
         x = "(log) Non-temptation Consumption [R$]",
         y = "(log) Temptation Consumption [R$]")
    
  
  plot4 <-  ggplot(temptation_cross, aes(x = log(`non-temptation_pc`), y = log(temptation_pc))) +
    geom_point(size = 0.5) +
    geom_smooth(model = 'lm', formula = y ~ x +1, se = FALSE, col = "red") +
    labs(title = "Temptation vs Non-temptation consumption",
         subtitle = "POF 2017 - 2018 (per-capita)",
         x = "(log) Non-temptation Consumption [R$]",
         y = "(log) Temptation Consumption [R$]") +
    scale_x_continuous(breaks = seq(2,14,2)) +
    scale_y_continuous(breaks = seq(2,14,2))
  
  p <- plot_grid(plot2,plot4)
  
  if(save){ 
    
      ggsave(glue("Temptation - {ex_name} - fraction of consumption, scatter.png"),
             plot = plot1,
             width = 5, height = 4.5)
      ggsave(glue("Temptation - {ex_name} - fraction of consumption, smooth.png"),
             plot = plot2,
             width = 5, height = 4.5)
      ggsave(glue("Temptation - {ex_name} - Temptation by non-temptation, smooth.png"),
             plot = plot3,
             width = 5, height = 4.5)
      ggsave(glue("Temptation - {ex_name} - scatter and lm.png"),
             plot = plot4,
             width = 5, height = 4.5)
      
      ggsave(glue("Temptation - {ex_name} - scatter and lm sbs fraction smooth.png"), 
             plot = p, 
             width = 12, height = 5)
  }
}

#-----------------------------------------------------------------------
#                     DEFINING EXERCISE SETTINGS AND RUNNING
#-----------------------------------------------------------------------

trim_outliers <- FALSE
lower_trim <- 0       # monthly per capita
upper_trim <- 20000     # monthly per capita

remove_durables <- TRUE

# Running exercise:

sigma <-  regress_and_plot(consumption_total = consumption_total,
                           morador_count = morador_count,
                           ex_name = "No Leisure, 20k Trim",
                           trim_outliers = trim_outliers,
                           lower_trim = lower_trim,
                           upper_trim = upper_trim,
                           remove_durables = remove_durables,
                           save = TRUE)

#-----------------------------------------------------------------------
#                     ADDITIONAL PLOTS
#-----------------------------------------------------------------------

setwd(export_path)

ggplot(data = temptation_filtered) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac))

ggplot(data = temptation_filtered %>% filter(total_monthly_pc<5000)) +
  geom_histogram(aes(x = total_monthly_pc), breaks = seq(0,5000,100))

ggplot(temptation_filtered) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
  labs(title = "Temptation by Consumption level",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  coord_cartesian(xlim = c(0,3000), ylim = c(0, 0.12)) +
  ggsave(glue("Temptation - 96perc No Leisure - fraction of consumption, smooth.png"),
         width = 5, height = 4.5)

#-------------------------- Making a 96 percentile zoom --------------------

ggplot(temptation_cross) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
  labs(title = "Temptation by Consumption level",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  coord_cartesian(xlim = c(0,3000), ylim = c(0, 0.08)) +
  ggsave(glue("Temptation - 96perc, No Leisure, 20k Trim - fraction of consumption, smooth.png"),
         width = 5, height = 4.5)

#----------------------- Running Baseline regression -----------------------

temptation_filtered <- temptation_cross %>% 
  filter(temptation>0 , total_monthly_pc < upper_trim)

lm <- lm(log(temptation_pc) ~ log(`non-temptation_pc`) + 1, temptation_filtered)
summary(lm)

# Generating Stargazer Table
stargazer(lm, ci = T, style = "qje")

# Implied sigmas (under xi = 1)
gammas <- lm$coefficients

sigma_x_eq <- function(sigma_x, gamma1, gamma2, xi){
  error <- gamma2*sigma_x*log(gamma2*(sigma_x-1)/(gamma2*sigma_x - 1)*(1/xi)) + gamma1
  return(error^2)
}

optim <- optimize(interval = c(1.01,100),
                  f = sigma_x_eq,
                  gamma1 = gammas[1],
                  gamma2 = gammas[2],
                  xi = 2)

sigma_x <- optim$minimum
sigma_y <- gammas[2]*sigma_x
names(sigma_y) <- NULL

ggplot(temptation_filtered) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
  labs(title = "Temptation by Consumption level (excluding zeros)",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  coord_cartesian(xlim = c(0,20000), ylim = c(0, 0.12)) +
  ggsave(glue("Temptation - No-Zero Tempt - fraction of consumption, smooth.png"),
         width = 5, height = 4.5)

ggplot(temptation_filtered, aes(x = log(`non-temptation_pc`), y = log(temptation_pc))) +
  geom_point(size = 0.5) +
  geom_smooth(model = 'lm', formula = y ~ x +1, se = FALSE, col = "red") +
  labs(title = "Temptation vs Non-temptation consumption (excluding zeros)",
       x = "(log) Non-temptation Consumption [R$]",
       y = "(log) Temptation Consumption [R$]") +
  scale_x_continuous(breaks = seq(2,14,2)) +
  scale_y_continuous(breaks = seq(2,14,2)) +
  coord_cartesian(xlim = c(4,14), ylim = c(1, 11)) +
  ggsave(glue("Temptation - No-Zero Tempt - scatter and lm, smooth.png"),
         width = 6, height = 5.4)


#-------------- Checking if zero temptation has any trend ------------------

temptation_binary <- temptation_cross %>%
  mutate(has_tempt = temptation != 0)

lm_bin <- lm(has_tempt ~ total_monthly_pc + 1, temptation_binary)
summary(lm_bin)

ggplot(temptation_binary) +
  geom_smooth(aes(x = total_monthly_pc, y = as.numeric(has_tempt))) +
  coord_cartesian(xlim = c(0,3000), ylim = c(0, 1)) +
  labs(title = "Local Fraction of Sample with Positive Temptation",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Local Prob. of positive temptation consumption") +
ggsave(glue("Temptation - 96perc, No Leisure - local prob of positive tempt, smooth.png"),
       width = 5, height = 4.5)

# selecting < 3000 BRL/pc month sample
lm_bin <- lm(has_tempt ~ total_monthly_pc + 1, temptation_binary %>% filter(total_monthly_pc<3000))
summary(lm_bin)

# selecting 500 < x < 3000 BRL/pc month sample
lm_bin <- lm(has_tempt ~ total_monthly_pc + 1, temptation_binary %>% filter(total_monthly_pc<3000, total_monthly_pc>500))
summary(lm_bin)

# selecting 1000 < x < 3000 BRL/pc month sample
lm_bin <- lm(has_tempt ~ total_monthly_pc + 1, temptation_binary %>% filter(total_monthly_pc<3000, total_monthly_pc>1000))
summary(lm_bin)
