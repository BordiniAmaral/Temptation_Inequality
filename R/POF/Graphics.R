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

export_path <- "D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/POF/R Project/Exports"
setwd(export_path)

#-----------------------------------------------------------------------
#                     TIDYING EMPIRICAL TEMPTATION DF
#-----------------------------------------------------------------------

# General parameter for selected sample

remove_durables <- TRUE
trim_outliers <- TRUE
lower_trim <- 0       # monthly per capita
upper_trim <- 20000   # monthly per capita

# Actual tidying up

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

temptation_filtered <- temptation_cross %>%
  filter(total >0)

if(trim_outliers){
  temptation_filtered <- temptation_filtered %>%
    filter(total_monthly_pc > lower_trim, total_monthly_pc < upper_trim)
}

#-----------------------------------------------------------------------
#                     PLOTS
#-----------------------------------------------------------------------

setwd(export_path)

#-------------------------- Making a 96 percentile zoom --------------------

ggplot(temptation_cross) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
  labs(title = "Temptation by Consumption level",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  coord_cartesian(xlim = c(0,3000), ylim = c(0, 0.08)) +
  ggsave(glue("Temptation - 96perc - fraction of consumption, smooth.png"),
         width = 5, height = 4.5)

#----------------------- Baseline Graphs ----------------------------------

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
ggsave(glue("Temptation - 96perc - local prob of positive tempt, smooth.png"),
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
