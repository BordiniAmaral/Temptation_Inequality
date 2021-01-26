#-----------------------------------------------------------------------
# This file is dedicated to exploring the stoney-geary formulation in an
# environment with temptation.
#
# For several x0 guesses, run regressions on intratemporal optimality
# conditions, obtain implied sigmas and plots simulated curves
#
# [REQUIRES]
# consumption_total dataframe obtained from tidying_data.R
# morador_count dataframe obtained from tidying_data.R
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
library(dtplyr)

#-----------------------------------------------------------------------
#                     DEFINING EXPORTING PATH
#-----------------------------------------------------------------------

# WARNING #
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

fwrite(temptation_filtered, "temptation_filtered.csv", sep = ";")

#-----------------------------------------------------------------------
#      Plotting curves from gammas and x0 (for arbitrary tests)
#-----------------------------------------------------------------------

y_sim_gamma <- function(x0,gamma1,gamma2,x){

  # Filtering relevant x > x0
  relevant <- x > x0
  x <- x[relevant]
  y <- exp(gamma1 + gamma2*log(x-x0))

  results <- data.frame(x = x,
                        y = y)

  results <- results %>%
    mutate(total = x + y,
           tempt_frac = y / total,
           x0 = x0)

  return(results)
} # REMEMBER TO CHECK IF FUNCTIONAL FORM IS UP TO DATE

x_points <- seq(0,240000,100) # DEFINE X-AXIS RANGE

settings <- list(x0 = 0,
                 gamma1 = 1,
                 gamma2 = 1) # INSERT VALUES HERE

teste <- pmap(settings, y_sim_gamma, x = x_points) %>% bind_rows()

ggplot() +
  geom_line(data = teste, aes(x = total/12, y = tempt_frac, colour = factor(x0))) +
  geom_smooth(data = temptation_filtered, aes(x = total_monthly_pc, y = tempt_frac))

ggplot() +
  geom_point(data = temptation_filtered, aes(x = total_monthly_pc, y = tempt_frac), size = 0.7) +
  geom_line(data = teste, aes(x = total/12, y = tempt_frac), colour = "red", size = 1)


#-----------------------------------------------------------------------
#   Running regressions (Cases A,C), inserting GMM result (Case B)
#-----------------------------------------------------------------------

# Function for obtaining implied sigmas from gammas, under arbitrary xi
sigma_x_eq <- function(sigma_x, gamma1, gamma2, xi){
  error <- gamma2*sigma_x*log(gamma2*(sigma_x-1)/(gamma2*sigma_x - 1)*(1/xi)) + gamma1
  return(error^2)
} # REMEMBER TO CHECK IF FUNCTIONAL FORM IS UP TO DATE

#--------------------------- Case A ------------------------------------

# Define x0 (annual) used: 0 for case A; the same as case B for case C
x0 <- 75*12
upper_trim <- 20000 # monthly per capita

temptation_filtered_regA <- temptation_cross %>%
  filter(temptation > 0 , total_monthly_pc < upper_trim, total_monthly_pc*12 > x0)

lm_A <- lm(log(temptation_pc) ~ log(`non-temptation_pc` - x0) + 1, temptation_filtered_regA)
summary(lm_A)

# Generating Stargazer Table
stargazer(lm_A, ci = T, style = "qje") # SEE IF IT IS POSSIBLE TO IMPROVE THIS OUTPUT

# Obtaining implied sigmas from gammas, under arbitrary xi (to use arbitrary sigma_x in the future)
gammas_A <- lm_A$coefficients
xi <- 1

optim <- optimize(interval = c(1.01,100),
                  f = sigma_x_eq,
                  gamma1 = gammas_A[1],
                  gamma2 = gammas_A[2],
                  xi = xi)

sigma_x_A <- optim$minimum
sigma_y_A <- gammas[2]*sigma_x
names(sigma_y_A) <- NULL

#--------------------------- Case B ------------------------------------

# Manually input these values from python GMM script

x0_B <- 0           # INSERT VALUE HERE
gammas_B <- c(1,2)  # INSERT VALUE HERE

# Obtaining implied sigmas from gammas, under arbitrary xi (to use arbitrary sigma_x in the future)
xi <- 1

optim <- optimize(interval = c(1.01,100),
                  f = sigma_x_eq,
                  gamma1 = gammas_A[1],
                  gamma2 = gammas_A[2],
                  xi = xi)

sigma_x_B <- optim$minimum
sigma_y_B <- gammas[2]*sigma_x
names(sigma_y_B) <- NULL

#--------------------------- Case C ------------------------------------


# Define x0 (annual) used: 0 for case A; the same as case B for case C
x0 <- 75*12
upper_trim <- 20000 # monthly per capita

temptation_filtered_regC <- temptation_cross %>%
  filter(temptation > 0 , total_monthly_pc < upper_trim, total_monthly_pc*12 > x0)

lm_C <- lm(log(temptation_pc) ~ log(`non-temptation_pc` - x0) + 1, temptation_filtered_regC)
summary(lm_C)

# Generating Stargazer Table
stargazer(lm_C, ci = T, style = "qje") # SEE IF IT IS POSSIBLE TO IMPROVE THIS OUTPUT

# Obtaining implied sigmas from gammas, under arbitrary xi (to use arbitrary sigma_x in the future)
gammas_C <- lm_C$coefficients
xi <- 1

optim <- optimize(interval = c(1.01,100),
                  f = sigma_x_eq,
                  gamma1 = gammas_C[1],
                  gamma2 = gammas_C[2],
                  xi = xi)

sigma_x_A <- optim$minimum
sigma_y_A <- gammas[2]*sigma_x
names(sigma_y_A) <- NULL


# To have a look, adapt which dataset you are selecting:

ggplot(temptation_filtered_reg) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
  labs(title = "Temptation by Consumption level (excluding zeros)",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  coord_cartesian(xlim = c(0,20000), ylim = c(0, 0.12))

ggplot(temptation_filtered_reg, aes(x = log(`non-temptation_pc` - x0), y = log(temptation_pc))) +
  geom_point(size = 0.5) +
  geom_smooth(model = 'lm', formula = y ~ x +1, se = FALSE, col = "red") +
  labs(title = "Temptation vs Non-temptation consumption (excluding zeros)",
       x = "(log) Non-temptation Consumption (above x0) [R$]",
       y = "(log) Temptation Consumption [R$]") +
  scale_x_continuous(breaks = seq(2,14,2)) +
  scale_y_continuous(breaks = seq(2,14,2)) +
  coord_cartesian(xlim = c(4,14), ylim = c(1, 11))

#-----------------------------------------------------------------------
#                   Comparing solutions A, B and C
#-----------------------------------------------------------------------

# Solution A: from regression with x0

settings <- list(x0 = 0,
                 gamma1 = gammas_A[1],
                 gamma2 = gammas_A[2])

testeA <- pmap(settings, y_sim_gamma, x = x_points) %>% bind_rows()

# Solution B: from GMM

settings <- list(x0 = x0_B,
                 gamma1 = gammas_B[1],
                 gamma2 = gammas_B[2])

testeB <- pmap(settings, y_sim_gamma, x = x_points) %>% bind_rows()

# Solution C: regression with x0 from solution B

settings <- list(x0 = x0_B,
                 gamma1 = gammas_C[1],
                 gamma2 = gammas_C[2])

testeC <- pmap(settings, y_sim_gamma, x = x_points) %>% bind_rows()

# Comparative graph for tempt frac
ggplot() +
  geom_line(data = testeA, aes(x = total/12, y = tempt_frac, color = "(A)"), linetype = "dashed", size = 1.2) +
  geom_line(data = testeB, aes(x = total/12, y = tempt_frac, color = "(B)"), linetype = "dashed", size = 1.2) +
  geom_line(data = testeC, aes(x = total/12, y = tempt_frac, color = "(C)"), linetype = "dashed", size = 1.2) +
  geom_smooth(data = temptation_filtered, aes(x = total_monthly_pc, y = tempt_frac, color = "Data")) +
  # geom_vline(xintercept = gridm_v, linetype = "dashed", linesize = 0.5) +
  coord_cartesian(xlim = c(0,3000), ylim = c(0, 0.12)) +
  scale_color_manual(values = c(
    'Data' = 'blue',
    '(A)' = 'black',
    '(B)' = 'red',
    '(C)' = 'yellow')) +
  labs(color = '') +
  labs(title = "Calibrated Curves",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  theme_minimal() +
  theme(legend.text = element_text(size = 16),
        legend.key.size = unit(1.5,"cm"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  ggsave(glue("Calibrations_ABC_Tempt_frac96.png"), width = 10, height = 6.5)

# Comparative graph for log(tempt) vs log(non-tempt)

ggplot() +
  geom_point(data = temptation_filtered, aes(x = log(`non-temptation_pc`), y = log(temptation_pc)), size = 0.7, colour = "grey51") +
  geom_smooth(data = temptation_filtered, aes(x = log(`non-temptation_pc`), y = log(temptation_pc), colour = "(A)"), size = 1.2, model = 'lm', formula = y ~ x +1, se = FALSE) +
  geom_line(data = testeB, aes(x = log(x), y = log(y), colour = "(B)"), size = 1.2)+
  geom_line(data = testeC, aes(x = log(x), y = log(y), colour = "(C)"), size = 1.2)+
  coord_cartesian(xlim = c(4,14), ylim = c(2, 10)) +
  labs(title = "Temptation by Consumption level - Calibrations (A), (B) and (C)",
       x = "(log) Non-temptation Consumption [R$]",
       y = "(log) Temptation Consumption [R$]")+
  scale_color_manual(values = c(
    '(A)' = 'black',
    '(B)' = 'red',
    '(C)' = 'yellow')) +
  labs(color = '') +
  theme_minimal() +
  theme(legend.text = element_text(size = 16),
        legend.key.size = unit(1.5,"cm"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  scale_x_continuous(breaks = seq(2,14,2)) +
  scale_y_continuous(breaks = seq(2,14,2)) +
  ggsave(glue("Temptation - No-Zero Tempt - scatter and line, three alternatives.png"), width = 10, height = 6.5)

