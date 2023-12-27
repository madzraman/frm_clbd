# Assume that HOMOGENEOUS data is generated from sas in one batch of 250 datasets and another batch of 750 datasets for each value of Beta_WS.

library(lme4)
library(survival)
library(tidyverse)
library(haven)

num_datasets_each <- 250
num_runs <- 4
num_subjects <- 200
num_timepoints <- 10
num_param_combos <- 9 # number of beta_BS and cov combos

n <- integer(num_datasets_each*num_subjects*num_timepoints)

for (i in 1:num_datasets_each){
    z <- rep(i, num_subjects*num_timepoints) + (750) # because I think we had a 250 set and a 750 set 
    n[((i-1)*(num_subjects*num_timepoints) + 1):(i*(num_subjects*num_timepoints))] <- z
}

homo_out_250d_n04$ndat <- rep(n, num_param_combos)

simdat <- rbind(homo_out_750d_n04, homo_out_250d_n04)
