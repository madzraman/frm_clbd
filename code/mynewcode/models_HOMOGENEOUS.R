# Assume data is already processed and combined into the one large dataset.
# Specifically, the input of this file is the output of data_processing_HOMOGENEOUS.R

library(lme4)
library(survival)
library(tidyverse)
library(haven)

# Define constants based on SAS data table and reduce df size (reduce number of columns)

ndats <- 1000 #simdat$ndats[1], 1000
npergrp <- simdat$npergrp[1] # then 2 groups
nsubj <- simdat$nsubj[1] # so 200 subjects total, but we aren't really using/dividing the groups here
ssbetaWS <- simdat$ssbeta[1] #change back to ssbetaWS
all_ssbetaBS <- c(simdat$ssbetaBS1[1], simdat$ssbetaBS2[1], simdat$ssbetaBS3[1]) # k
ntimepoints <- simdat$ntime[1]
icc <- simdat$icc[1]
all_covs <- c(simdat$cov1[1], simdat$cov2[1], simdat$cov3[1]) # j
varerr <- pi^2 / 3  # note that atan(1) = pi/4 so atan(1) * 4 = pi. varerr = (atan(1) * 4)**2 / 3
varsub <- varerr * (icc / (1-icc)) # subject variance given the icc and error variance (sigma squared v given icc and sigma squared e).
sdsub <- sqrt(varsub)

# Keep only the necessary columns for running the models
simdat_new <- simdat[,-c(1:16, 20)] 


generate_blank_results_matrices <- function(){
    all_estimates_mat <- matrix(nrow = ndats, ncol = 9)#length(all_covs)*length(all_ssbetaBS))
    cov_rej_mat <- matrix(0, nrow = 2, ncol = 9)#length(all_covs)*length(all_ssbetaBS))
    return(list(all_estimates_mat, cov_rej_mat))
}

poss_model_names <- c("randint", "randintBSWS", "clogit", "randintWS")

get_model_results <- function(dataset_number, condition_number, model_name, model_output, all_estimates_matrix, cov_rej_matrix){
    if (model_name == "randint") {
        est <- summary(model_output)$coefficients[2,1]
        stderr <- summary(model_output)$coefficients[2,2]
    }
    else if (model_name == "randintBSWS"){
        est <- summary(model_output)$coefficients[3,1]
        stderr <- summary(model_output)$coefficients[3,2]
    }
    else if (model_name == "clogit"){
        est <- summary(model_output)$coefficients[1,1]
        stderr <- summary(model_output)$coefficients[1,3]
    }
    else if (model_name == "randintWS"){
        est <- summary(model_output)$coefficients[2,1]
        stderr <- summary(model_output)$coefficients[2,2]
    }
    m.upper <- est + qnorm(0.975)*stderr
    m.lower <- est - qnorm(0.975)*stderr
    covers_inrange <- ifelse((ssbetaWS > m.lower) & (ssbetaWS < m.upper), 1, 0)
    reject_signif <- ifelse((0 < m1.lower) | (0 > m1.upper), 1, 0)
    all_estimates_matrix[dataset_number, condition_number] <- est
    cov_rej_matrix[1, condition_number] <- cov_rej_matrix[1, condition_number] + covers_inrange # coverage
    cov_rej_matrix[2, condition_number] <- cov_rej_matrix[2, condition_number] + reject_signif # rejection
    
}

m1.all_estimates_mat <- generate_blank_results_matrices()[[1]]
m1.cov_rej_mat <- generate_blank_results_matrices()[[2]]
m2.all_estimates_mat <- generate_blank_results_matrices()[[1]]
m2.cov_rej_mat <- generate_blank_results_matrices()[[2]]
m3.all_estimates_mat <- generate_blank_results_matrices()[[1]]
m3.cov_rej_mat <- generate_blank_results_matrices()[[2]]
m4.all_estimates_mat <- generate_blank_results_matrices()[[1]]
m4.cov_rej_mat <- generate_blank_results_matrices()[[2]]


for (dd in 1:ndats){ # for each dataset
    if (dd %% 100 == 0){ # counter 
        print(dd)
    }
    this_specific_df <- simdat_new[simdat_new$ndat == dd,]
    condition_num <- 0
    # j, k, ndat are indicators of the loop for which ssbetaBS and which cov and which dataset we're on. 
    # we repeat the 9 combinations for 10 datasets each. 
    for (jj in 1:length(all_covs)){ 
        for (kk in 1:length(all_ssbetaBS)){
            condition_num <- condition_num + 1 # 1 through 9
            this_df <- this_specific_df[(this_specific_df$j == jj) & (this_specific_df$k == kk),]
            this_cov <- all_covs[jj]
            this_ssbetaBS <- all_ssbetaBS[kk]
            
            # Homogeneous exposure models
            m1.out <- glmer(y ~ grp + (1|id), family=binomial, data=this_df) 
            m2.out <- glmer(y ~ grpmean + grpdev + (1|id), family=binomial, data=this_df) # grpdev = grp - grpmean
            m3.out <- clogit(y ~ grp + strata(id), data=this_df) # conditional logistic regression model
            m4.out <- glmer(y ~ grpdev + (1|id), family=binomial, data=this_df) # adaptive centering approach, WS effect only
            
            # M1: Random Intercept model assuming BS=WS effect, no BS/WS decomposition
            get_model_results(dd, condition_num, "randint", m1.out, m1.all_estimates_mat, m1.cov_rej_mat)
 
            # M2: Random Intercept model with BS/WS decomposition
            get_model_results(dd, condition_num, "randintBSWS", m2.out, m2.all_estimates_mat, m2.cov_rej_mat)
            
            # M3: Conditional Logistic Regression
            get_model_results(dd, condition_num, "clogit", m3.out, m3.all_estimates_mat, m3.cov_rej_mat)

            # M4: Random Intercept model with Adaptive Centering approach (WS deviation effect only)
            get_model_results(dd, condition_num, "randintWS", m4.out, m4.all_estimates_mat, m4.cov_rej_mat)
        }
    }
}



