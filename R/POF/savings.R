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

#-----------------------------------------------------------------------
#                     Stats by quantile group
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
