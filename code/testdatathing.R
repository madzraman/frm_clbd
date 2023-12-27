library(lme4)
library(survival)
library(tidyverse)
library(haven)


d0 <- data.frame(read_sas("/Users/madhuriraman/Downloads/dec_23_hetero_out_100d_n04_4.sas7bdat", NULL))
