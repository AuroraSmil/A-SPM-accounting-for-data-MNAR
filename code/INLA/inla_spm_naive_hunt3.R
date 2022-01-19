library(INLA)
library(brinla)
#INLA:::inla.dynload.workaround()
library(tidyverse)
library(data.table)
library(latex2exp)
library(tikzDevice)
library(cowplot)
library(gridExtra)
library(grid)
library(gridtext)
library(cowplot)
library(ggplot2)
library(glue)
library(dplyr)
library(DescTools)
library(scoringRules)
library(mgcv)


run_spm_H3 <- function(data_sim, c_mu = 0, c_sigma = 1, verbose = F, compute_all = FALSE){
  
  sigma2 = 0.001
  #Prepare the response variables:
  n <- nrow(data_sim)
  y_gaussian <- c(data_sim$bp_4_corr,rep(NA,n))
  m_binomial <- c(rep(NA,n),data_sim$missing)
  joint_response <- list(y_gaussian,m_binomial)
  
  #Must specify number of trials for binomial repsonse variable
  ntrials <- c(rep(NA,n),rep(1,n))
  
  linear_covariates <- data.frame(alpha_0 = c(rep(1,n),rep(NA,n)),
                                  beta_0 = c(rep(NA,n),rep(1,n)),
                                  y_SEX = c(data_sim$sex,rep(NA,n)),
                                  y_BP3 = c(data_sim$bp_3_corr,rep(NA,n)),
                                  y_AGE = c(data_sim$age_3,rep(NA,n)),
                                  y_BMI = c(data_sim$bmi_3,rep(NA,n)),
                                  m_SEX = c(rep(NA,n),data_sim$sex),
                                  m_BP3 = c(rep(NA,n),data_sim$bp_3_corr),
                                  m_AGE = c(rep(NA,n),data_sim$age_3),
                                  m_BMI = c(rep(NA,n),data_sim$bmi_3)
  )
  random_covariates <- data.frame(y_eps1 = c(1:n,rep(NA,n)),
                                  m_eps1 = c(rep(NA,n),1:n))
  
  joint_data <- c(linear_covariates,random_covariates)
  
  joint_data$Y <- joint_response
  
  formula_spm = Y ~ -1 + alpha_0 + y_BP3 + y_AGE + y_SEX + y_BMI +
    beta_0 + m_SEX + m_BP3 +  m_BMI +
    f(m_AGE,model="rw2",constr=T) +
    f(y_eps1,model="iid") +
    f(m_eps1,copy="y_eps1",fixed=F,param=c(c_mu,c_sigma)) #Copy epsilon1 into binomial dropout process with association parameter
  # Use the option param to set mean and precision on Gaussian prior for associaction parameter
  
  link= rep(NA, 2*(n ))
  link[which(is.na(y_gaussian[1:(n )]))] = 1
  link[(n) + which(is.na(m_binomial[((n)+1):(2*(n))]))] = 2
  length(link)
  start_time <- Sys.time()
  if(compute_all == FALSE ){
    joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                       Ntrials=ntrials, verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                       control.predictor = list(link = link))
  }else{
    joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                       Ntrials=ntrials, verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                       control.predictor = list(link = link),
                       control.compute=list(config=TRUE))
    
  }
  
  
  # Use control.family to fix sigma2. Specified through log-precision.
  
  end_time <- Sys.time()
  return(joint_inla)
}

run_naive_bp_H3 <- function(data_sim, verbose = F,compute_all = FALSE){
  n <- nrow(data_sim)
  sigma2 = 0.001 #Fixed to a neglectable value
  
  #Prepare the response variables:
  y_gaussian <- data_sim$bp_4_corr
  
  
  linear_covariates <- list(alpha_0 = rep(1,n),
                            y_SEX = data_sim$sex,
                            y_BP1 = data_sim$bp_3_corr,
                            y_AGE = data_sim$age_3,
                            y_BMI = data_sim$bmi_3,
                            y_eps1 = 1:n)
  naive_data_bp <- linear_covariates
  
  naive_data_bp$Y <- y_gaussian
  
  formula_naive_bp = Y ~ -1 + alpha_0 + y_BP1 + y_AGE + y_SEX + y_BMI +
    f(y_eps1,model="iid") 
  
  start_time <- Sys.time()
  if(compute_all == FALSE){
    joint_inla <- inla(formula_naive_bp, 
                       family="gaussian",
                       data=naive_data_bp, 
                       verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T)))
  }else{
    joint_inla <- inla(formula_naive_bp, 
                       family="gaussian",
                       data=naive_data_bp, 
                       verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T)),
                       control.compute=list(config=TRUE))
  }
  
  return(joint_inla)
}

run_naive_m_H3 <- function(data_sim, verbose = F, compute_all = FALSE){
  data<- data_sim
  n <- nrow(data)
  
  linear_covariates_n_m <- list(beta_0 = rep(1,n),
                                m_SEX = data$sex,
                                m_BP1 = data$bp_3_corr,
                                m_AGE = data$age_3,
                                m_BMI = data$bmi_3,
                                m_eps1 =1:n)
  
  naive_data_m <- linear_covariates_n_m
  y_naive_m <- data$missing
  
  naive_data_m$Y <- y_naive_m
  formula_naive_m = Y ~ -1 + beta_0  + m_SEX + m_BP1 + m_BMI + f(m_AGE,model="rw2",constr=T)
  ntrials <- 1
  
  if(compute_all == FALSE){
    inla_naive_m <- inla(formula = formula_naive_m, 
                         family = "binomial", 
                         Ntrials = ntrials, 
                         data = naive_data_m, 
                         verbose=verbose)
  }else{
    inla_naive_m <- inla(formula = formula_naive_m, 
                         family = "binomial", 
                         data = naive_data_m, 
                         Ntrials=ntrials, 
                         verbose=verbose,
                         control.compute=list(config=TRUE))
  }
  
  
  
  
  # Use control.family to fix sigma2. Specified through log-precision.
  
  return(inla_naive_m)
}

source("code/explore/explore_functions.R")
source("code/INLA/compare_spm_naive_estimates_hunt2_cohort.R")

# Set directory for storing the SPM and naive model
dir <- "/home/aurorach/data_HUNT_aurora/master/"
# Read data ---- 
data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")
data_small <- as.data.table(data)#[1:1000])
n <- nrow(data_small)
data <- data_small

sigma2 = 0.001 #Fixed to a neglectable value

# Run SPM ---- 

inla_spm <- run_spm_H3(data, verbose = T, compute_all = FALSE)


saveRDS(inla_spm, file = glue::glue("{dir}spm_inla_HUNT34.RData"))
summary(inla_spm)

# Run naive -----
inla_naive_bp <- run_naive_bp_H3(data, verbose = T, compute_all = FALSE)
saveRDS(inla_naive_bp, file = glue::glue("{dir}inla_naive_bp_HUNT34.RData"))

inla_naive_m <- run_naive_m_H3(data, verbose = T, compute_all = FALSE)
saveRDS(inla_naive_m, file = glue::glue("{dir}inla_naive_m_HUNT34.RData"))

summary(inla_naive_bp)
summary(inla_naive_m)

# Compare SPM and naive model estimates ----
# Can be run without the above code, just remember to load the packages.

compare_spm_estimates(dir, cohort_number = 3)


# inla_naive_bp <- readRDS( file = glue::glue("{dir}inla_naive_bp_HUNT34.RData"))
# inla_naive_bp$summary.fixed
# bri.hyperpar.summary(inla_naive_bp)
# 
# inla_naive_m <- readRDS( file = glue::glue("{dir}inla_naive_m_HUNT34.RData"))
# inla_naive_m$summary.fixed
# bri.hyperpar.summary(inla_naive_m)
# 
# h3_fixed_naive <- c(inla_naive_bp$summary.fixed$mean, 
#                     inla_naive_m$summary.fixed$mean,
#                     bri.hyperpar.summary(inla_naive_bp)[1], 
#                     bri.hyperpar.summary(inla_naive_m)[1] )
# inla_naive_bp <- readRDS( file = glue::glue("{dir}inla_naive_bp_HUNT34.RData"))
# inla_naive_bp$summary.fixed
# bri.hyperpar.summary(inla_naive_bp)
# 
# inla_naive_bp_h2 <- readRDS( file = glue::glue("{dir}inla_naive_bp.RData"))
# inla_naive_m_h2 <- readRDS( file = glue::glue("{dir}inla_naive_m.RData"))
# 
# h2_fixed_naive <- c(inla_naive_bp_h2$summary.fixed$mean, 
#                     inla_naive_m_h2$summary.fixed$mean,
#                     bri.hyperpar.summary(inla_naive_bp_h2)[1], 
#                     bri.hyperpar.summary(inla_naive_m_h2)[1] )
# 
# h2_fixed_naive - h3_fixed_naive
# 
# inla_spm_h2 <- readRDS(glue::glue("{dir}spm_inla_HUNT23.RData"))
# inla_spm_h3 <- readRDS(glue::glue("{dir}spm_inla_HUNT34.RData"))
# 
# inla_spm_h2$summary.fixed$mean - inla_spm_h3$summary.fixed$mean
# 
# bri.hyperpar.summary(inla_spm_h2)[1:3]- bri.hyperpar.summary(inla_spm_h3)[1:3]
