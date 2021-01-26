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

# Cleaning undesired objects
rm(list = setdiff(ls(), c("consumption_total","morador_count")))

export_path <- "C:/Backup Acer/EESP-Mestrado/Dissertação/POF/R Project/Exports"
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

#########################################################################
##                     2-steps section (so far, unsuccessful)          ##
#########################################################################
#-----------------------------------------------------------------------
#                   Regressing for several x0
#-----------------------------------------------------------------------

obtain_gammas_sg <- function(x0, temptation_filtered){
  
  cat("------------------------------\n")
  cat(glue("Regressing for x0 = {x0} \n\n"))
  temptation_sg <- temptation_filtered %>% 
    mutate(x_x0 = `non-temptation_pc` - x0) %>%
    filter(x_x0 > 0)
  
  lm <- lm(log(temptation_pc) ~ log(x_x0) + 1, temptation_sg)
  gamma <- lm$coefficients
  names(gamma) <- NULL
  n <- nrow(temptation_sg)
  
  results <- c(gamma[1],gamma[2],x0, n)
  return(results)
}

# Remember: all income is in an yearly basis
# Using for this example R$0 (no sg) to R$ 500 per capita monthly
x0 <- seq(0,500*12,10*12)

gammas <- map(x0, obtain_gammas_sg, temptation_filtered) %>%
  unlist() %>% 
  matrix(ncol = 4, byrow = TRUE) %>%
  data.frame() %>%
  rename(gamma1 = X1, gamma2 = X2, x0 = X3, observations = X4)

#-----------------------------------------------------------------------
#                   For each regression, testing several xi
#-----------------------------------------------------------------------

sigma_x_eq <- function(sigma_x, gamma1, gamma2, xi){
  error <- gamma2*sigma_x*log(gamma2*(sigma_x-1)/(gamma2*sigma_x - 1)*(1/xi)) + gamma1
  return(error^2)
}

obtain_sigmas_sg <- function(xi, gamma1, gamma2, x0, obs){
  
  cat("------------------------------\n")
  cat(glue(" Currently on x0 = {x0} and xi = {xi} \n\n"))
  
  optim <- optimize(interval = c(1.01,100),
                    f = sigma_x_eq,
                    gamma1 = gamma1,
                    gamma2 = gamma2,
                    xi = xi)
  
  sigma_x <- optim$minimum
  sigma_y <- gamma2*sigma_x
  names(sigma_y) <- NULL
  
  results <- c(sigma_x, sigma_y, xi, x0, gamma1, gamma2, obs)
  return(results)
}

xi <- seq(0.1,20,0.1)

settings <- list(xi = rep(xi, length(gammas$x0)),
                 gamma1 = rep(gammas$gamma1, each = length(xi)),
                 gamma2 = rep(gammas$gamma2, each = length(xi)),
                 x0 = rep(gammas$x0, each = length(xi)),
                 obs = rep(gammas$observations, each = length(xi)))

sigmas <- pmap(settings,obtain_sigmas_sg) %>%
  unlist() %>% 
  matrix(ncol = 7, byrow = TRUE) %>%
  data.frame() %>%
  rename(sigma_x = X1, sigma_y = X2, xi = X3, x0 = X4, gamma1 = X5, gamma2 = X6, obs = X7)

#-----------------------------------------------------------------------
#            Plotting simulated curves (need only gammas and x0)
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
}

x_points <- seq(0,240000,100)

settings <- list(x0 = gammas$x0,
                 gamma1 = gammas$gamma1,
                 gamma2 = gammas$gamma2)

teste <- pmap(settings, y_sim_gamma, x = x_points) %>% bind_rows()

ggplot() +
  geom_line(data = teste, aes(x = total/12, y = tempt_frac, colour = factor(x0))) +
  geom_smooth(data = temptation_filtered, aes(x = total_monthly_pc, y = tempt_frac))

ggplot() +
  geom_point(data = temptation_filtered, aes(x = total_monthly_pc, y = tempt_frac), size = 0.7) +
  geom_line(data = teste, aes(x = total/12, y = tempt_frac), colour = "red", size = 1)

#########################################################################
##                     Full-Space section                              ##
#########################################################################
#-----------------------------------------------------------------------
#                    Running MSE, in stages of approximation
#-----------------------------------------------------------------------

simulate_xy <- function(x0,sigma_x,sigma_y,xi,x,y,frac, pb){
  
  pb$tick()$print()
  
  # Filtering relevant x > x0
  relevant <- x > x0
  x <- x[relevant]
  y <- y[relevant]
  frac <- frac[relevant]
  y_sim <- exp((sigma_y/sigma_x)*log(x-x0) - sigma_y*log((sigma_x - 1)/(sigma_y-1)*(sigma_y/sigma_x)*(1/xi)))
  
  tot_deciles <- quantile(x+y, probs = seq(0.05,1,0.05))
  tot_values <- data.frame(x = x,
                           y = y,
                           frac = frac) %>%
    rowwise() %>%  
    lazy_dt() %>% 
    mutate(tot_decile = sum(tot_deciles < (x+y))) %>% 
    group_by(tot_decile) %>% 
    summarise(mean_frac = sum(y) / (sum(x)+sum(y)), .groups = "drop") %>% 
    as.data.frame()
  
  sim_values <- data.frame(x = x,
                           y = y,
                           frac = frac,
                           y_sim = y_sim) %>% 
    rowwise() %>% 
    lazy_dt() %>% 
    mutate(tot_decile = sum(tot_deciles < (x+y))) %>%
    group_by(tot_decile) %>% 
    summarise(sim_frac = sum(y_sim) / (sum(x)+sum(y_sim)), .groups = "drop") %>% 
    as.data.frame()
    
  MSE <- tot_values %>% 
    left_join(sim_values, by = "tot_decile") %>% 
    mutate(dif_sq = (sim_frac - mean_frac)^2) %>%
    summarise(MSE = sum(dif_sq)/n())
  
  MSE <- data.frame(x0 = x0,
                    sigma_x = sigma_x,
                    sigma_y = sigma_y,
                    xi = xi,
                    MSE = MSE[1])
  
  return(MSE)
}

run_simulations <- function(mapping_grid, data){
  
  # create the progress bar with a dplyr function. 
  pb <- progress_estimated(nrow(expand.grid(mapping_grid)))
  
  cat("Simulating... \n")
  simulations <- purrr::pmap(expand.grid(mapping_grid), simulate_xy,
                              x = data[["ntempt_pc"]],
                              y = data[["tempt_pc"]],
                              frac = data[["frac"]], pb) %>% bind_rows()
  return(simulations)
}

y_sim_par <- function(x, x0,sigma_x,sigma_y,xi){
  
  # Filtering relevant x > x0
  relevant <- x > x0
  x <- x[relevant]
  y <- exp((sigma_y/sigma_x)*log(x-x0) - sigma_y*log((sigma_x - 1)/(sigma_y-1)*(sigma_y/sigma_x)*(1/xi)))
  
  results <- data.frame(x = x,
                        y = y)
  
  results <- results %>% 
    mutate(total = x + y,
           tempt_frac = y / total,
           x0 = x0)
  
  return(results)
}


# Starting with a wide grid, small interval, to evaluate:

mapping_grid <- list(x0 = seq(0*12, 400*12, 25*12),
                     sigma_x = seq(2, 6, 0.1),
                     sigma_y = seq(1.4, 5, 0.1),
                     xi = seq(1, 3.5, 0.5))


data <- list(total_pc = temptation_filtered$total_monthly_pc *12,
             tempt_pc = temptation_filtered$temptation_pc,
             ntempt_pc = temptation_filtered$`non-temptation_pc`,
             frac = temptation_filtered$tempt_frac)

simulations <- run_simulations(mapping_grid, data)


save_best2 <- simulations %>% filter(MSE < quantile(simulations$MSE,probs = 0.1))
  
# Plotting some to see how they look
x_points <- seq(0,240000,100)

testeB2 <- y_sim_par(x = x_points,
                   x0 = 75*12,
                   sigma_x = 6.04,
                   sigma_y = 5.35,
                   xi = 0.78)

# gridm_v <- c(0, 2400, 4800, 7200, 9600, 12000, 16800, 21600, 26400, 31200, 49200, 67200, 85200, 103200, 121200, 181200) / 12
# gridm_v <- perc10


# Plotting only one solution
ggplot() +
  geom_line(data = testeC, aes(x = total/12, y = tempt_frac, color = "Alt. B"), linetype = "dashed", size = 0.8) +
  geom_smooth(data = temptation_filtered, aes(x = total_monthly_pc, y = tempt_frac, color = "Data")) +
  # geom_vline(xintercept = gridm_v, linetype = "dashed", linesize = 0.5) +
  coord_cartesian(xlim = c(0,3000), ylim = c(0, 0.12)) +
  labs(color = '') +
  labs(title = "Calibrated Curves",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation")

# Plotting log x log consumption fit
ggplot() +
  geom_point(data = temptation_filtered, aes(x = log(`non-temptation_pc`), y = log(temptation_pc)), size = 0.7) +
  geom_line(data = testeC, aes(x = log(x), y = log(y)), colour = "red", size = 1)+
  coord_cartesian(xlim = c(4,14), ylim = c(2, 10)) +
  labs(title = "Temptation by Consumption level (SMM, full sample)",
       x = "(log) Non-temptation Consumption [R$]",
       y = "(log) Temptation Consumption [R$]")+
  scale_x_continuous(breaks = seq(2,14,2)) +
  scale_y_continuous(breaks = seq(2,14,2)) # +
  ggsave(glue("Temptation - No-Zero Tempt - scatter and line, SMM - full sample.png"), width = 6, height = 5.4)

#----------------------- Running Alternative C regression -----------------------

x0 <- 75*12

temptation_filtered_C2 <- temptation_cross %>% 
  filter(temptation > 0 , total_monthly_pc < upper_trim, total_monthly_pc*12 > x0)

lm <- lm(log(temptation_pc) ~ log(`non-temptation_pc` - x0) + 1, temptation_filtered_C2)
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


# Obtaining gammas from sigmas and xi (To study GMM results)
 sigx = 3.58
 sigy = 2.16
 xi = 2
 gamma1_inv = -sigy*log(((sigx-1)/sigx)*(sigy/(sigy-1))*(1/xi))
 gamma2_inv = sigy/sigx
    
 xi_change = 0.452
 optim <- optimize(interval = c(1.01,1000),
                   f = sigma_x_eq,
                   gamma1 = gamma1_inv, # or gammas[1]
                   gamma2 = gamma2_inv, # or gammas[2]
                   xi = xi_change)
 
 sigma_x <- optim$minimum
 sigma_y <- gamma2_inv*sigma_x # or gammas[2]
 gamma1_check <- -sigma_y*log(((sigma_x-1)/sigma_x)*(sigma_y/(sigma_y-1))*(1/xi_change))
 gamma2_check <-  sigma_y/sigma_x
 names(sigma_y) <- NULL

ggplot(temptation_filtered_C2) +
  geom_smooth(aes(x = total_monthly_pc, y = tempt_frac)) +
  labs(title = "Temptation by Consumption level (excluding zeros)",
       subtitle = "POF 2017 - 2018 (per-capita)",
       x = "Monthly expenditure per capita [R$]",
       y = "Fraction spent on temptation") +
  coord_cartesian(xlim = c(0,20000), ylim = c(0, 0.12))

ggplot(temptation_filtered_C2, aes(x = log(`non-temptation_pc` - x0), y = log(temptation_pc))) +
  geom_point(size = 0.5) +
  geom_smooth(model = 'lm', formula = y ~ x +1, se = FALSE, col = "red") +
  labs(title = "Temptation vs Non-temptation consumption (excluding zeros)",
       x = "(log) Non-temptation Consumption (above x0) [R$]",
       y = "(log) Temptation Consumption [R$]") +
  scale_x_continuous(breaks = seq(2,14,2)) +
  scale_y_continuous(breaks = seq(2,14,2)) +
  coord_cartesian(xlim = c(4,14), ylim = c(1, 11))

#----------------------- Comparing solutions A, B and C -----------------------

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

#-------------------- Rascunhera ------------------------


ggplot() +
  geom_point(data = temptation_filtered %>% filter(total_monthly_pc<5000), aes(x = total_monthly_pc, y = tempt_frac), size = 0.7) +
  geom_smooth(data = temptation_filtered %>% filter(total_monthly_pc<5000), aes(x = total_monthly_pc, y = tempt_frac))

ggplot() +
  geom_smooth(data = temptation_filtered %>% filter(total_monthly_pc<5000), aes(x = total_monthly_pc, y = tempt_frac))

ggplot() +
  geom_point(data = temptation_filtered, aes(x = log(`non-temptation_pc`), y = log(temptation_pc)), size = 0.5) +
  geom_line(data = teste, aes(x = log(x), y = log(y), colour = factor(x0)))

ggplot() +
  geom_smooth(data = temptation_filtered, aes(x = log(`non-temptation_pc`), y = log(temptation_pc)), size = 0.5) + 
  geom_line(data = teste, aes(x = log(x), y = log(y), colour = factor(x0)))

ggplot() +
  geom_line(data = teste %>% filter(total/12<5000), aes(x = total/12, y = tempt_frac, colour = factor(x0))) +
  geom_smooth(data = temptation_filtered %>% filter(total_monthly_pc<5000), aes(x = total_monthly_pc, y = tempt_frac))
