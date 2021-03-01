#-----------------------------------------------------------------------
# This file is dedicated to analyzing the savings rate based upon income
# and expenditures from POF 2018
#
# [REQUIRES]
# savings_by_residual object from tidying_data.R file
#
# File by: Felipe Bordini do Amaral
# Institution: FGV - EESP
# Project: Master Dissertation
#
#-----------------------------------------------------------------------
#                     STARTING NECESSARY PACKAGES
#-----------------------------------------------------------------------

library(tidyverse)
library(glue)
library(data.table)

#-----------------------------------------------------------------------
#                     Importing csv files from simulation [Must]
#-----------------------------------------------------------------------

# Insert simulation ID and import savings rate by quantile studied
sim = "xxx"
savings_nt <- read.csv(glue("D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/Temptation/Exports/Python (Post Jan21)/sim_{sim}_savings_nt.csv"), row.names = 1, header = TRUE, sep = ";")
savings_eq <- read.csv(glue("D:/Google Drive/Working Cloud/EESP-Mestrado/Dissertação/Temptation/Exports/Python (Post Jan21)/sim_{sim}_savings_eq.csv"), row.names = 1, header = TRUE, sep = ";")

#Tidying table
savings_nt <- savings_nt %>%
  bind_cols(as_tibble_col(rownames(savings_nt), column_name = "quant_below")) %>%
  pivot_longer(!quant_below, names_to = "age", values_to = "s_rate") %>%
  mutate(age = as.integer(str_sub(age,2)))

savings_eq <- savings_eq %>%
  bind_cols(as_tibble_col(rownames(savings_nt), column_name = "quant_below")) %>%
  pivot_longer(!quant_below, names_to = "age", values_to = "s_rate") %>%
  mutate(age = as.integer(str_sub(age,2)))

#-----------------------------------------------------------------------
#                     Preparing some data [Must]
#-----------------------------------------------------------------------

# Getting household age reference, adding to savings by residual
hh_age <- morador_aggr %>%
  filter(V0306 == 1) %>%
  select(hh_id,V0403) %>%
  rename(age = V0403)

savings_by_residual <- savings_by_residual %>%
  left_join(hh_age, by = c("id"="hh_id"))

#-----------------------------------------------------------------------
#                     Stats by quantile group [Optional]
#-----------------------------------------------------------------------

# Identifiying income quantile levels
income_quants <- savings_by_residual$income_pc %>%
  quantile(probs = c(seq(0,0.9,0.1),seq(0.92,1,0.02)), type = 8)

# Fetching names if necessary later
tops <- names(income_quants)

# Adding quantile info to table
savings_with_quantile <- savings_by_residual %>%
  rowwise() %>%
  mutate(quant = sum(income_quants < income_pc))

# Doing some summarised data on each quantile
savings_summarised <- savings_with_quantile %>%
  group_by(quant) %>%
  summarise(mean = sum(gross_savings_pc)/sum(income_pc),
            stdv = sd(savings_rate),
            wt_stdv = sqrt(sum(income_pc/sum(income_pc)*(savings_rate-mean)^2)))

#-----------------------------------------------------------------------
#                     Stats by age [Optional]
#-----------------------------------------------------------------------

# Creating quantiles
quant_by_five <- savings_by_residual %>%
  filter(age >= 25, age <=65) %>%
  mutate(rounded_five = floor(age/5)*5) %>%
  group_by(rounded_five) %>%
  summarise(q98 = quantile(income_pc, probs = 0.98, type = 8),
            q90 = quantile(income_pc, probs = 0.90, type = 8),
            q50 = quantile(income_pc, probs = 0.50, type = 8),
            q0 = 0) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quant")

# These below are not a very smart step, but it does the trick and I am tired
# to write them better

# Using more granular age measure
savings_by_age <- savings_by_residual %>%
  filter(age >= 25, age <=65) %>%
  mutate(rounded_five = floor(age/5)*5) %>%
  left_join(quant_by_five, by = "rounded_five") %>%
  group_by(id) %>%
  filter(income_pc > value) %>%
  filter(value == max(value)) %>%
  ungroup()

ggplot()+
  geom_smooth(data = savings_by_age, aes(x = age, y = savings_rate, color = quant), method = NULL)+
  geom_line(data = savings_eq, aes(x = age, y = s_rate, color = quant_below)) +
  geom_hline(yintercept = 0,color = 'black', linetype = "dashed") +
  coord_cartesian(xlim = c(25,65), ylim = c(-2, 1)) +
  scale_y_continuous(breaks = seq(-2,1,0.5), minor_breaks = seq(-3,1,0.1))

# Using ages grouped by five years
savings_by_five <- savings_by_residual %>%
  filter(age >= 25, age <=65) %>%
  mutate(rounded_five = floor(age/5)*5) %>%
  left_join(quant_by_five, by = "rounded_five") %>%
  group_by(id) %>%
  filter(income_pc > value) %>%
  filter(value == max(value)) %>%
  ungroup() %>%
  group_by(rounded_five, quant) %>%
  summarise(mean = sum(gross_savings_pc)/sum(income_pc),
            stdv = sd(savings_rate),
            wt_stdv = sqrt(sum(income_pc/sum(income_pc)*(savings_rate-mean)^2)))

ggplot(savings_by_five %>% filter(quant != "q0"))+
  geom_line(aes(x = rounded_five, y = mean, color = quant))


