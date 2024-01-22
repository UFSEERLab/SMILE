
# Testing the stochastic simulation

# Figures


stochastic.sims

library(tidyverse)

map_df(stochastic.sims, "LIZ")

stoch.s <- map_df(stochastic.sims, "Suceptibles") %>% 
  mutate(time = row_number()) %>% 
  pivot_longer(-time)

i<-1
my.s <- stochastic.sims[[i]]$Suceptibles
my.m <- stochastic.sims[[i]]$Immune
my.i <- stochastic.sims[[i]]$Infected
my.liz <- stochastic.sims[[i]]$LIZ
my.e <- stochastic.sims[[i]]$Environment

time <- 1:length(my.s)