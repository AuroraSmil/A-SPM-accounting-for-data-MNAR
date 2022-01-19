# generate 100 simulated dataset

library(tidyverse)
library(lubridate)
library(INLA)
library(data.table)
library(brinla)
library(ggplot2)
library(dplyr)
library(DescTools)
library(latex2exp)
library(tikzDevice)
library(scoringRules)
library(glue)
library(mgcv)
source("code/simulation_study/model_functions.R")
source("code/validation/validation_func.R")

# Create directory for storing results 
today <- today()
if(!dir.exists(glue("results/simulations/{today}"))){
  dir.create(glue("results/simulations/{today}"))
}
if(!dir.exists(glue("results/simulations/{today}/1"))){
  dir.create(glue("results/simulations/{today}/1"))
}

dir <- glue("results/simulations/{today}/1")


# Read models and data ---- 
inla_spm <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT23.RData")

dataH34<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")
#dataH23 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")
#number of simulations
n_sim = 100
# number of samples per simulation 
n_samples <- 300

#nH23 = nrow(dataH23)
nH34 = nrow(dataH34)
# Store spm parameters 

a_0  = inla_spm$summary.fixed[[1]][1]
a_bp = inla_spm$summary.fixed[[1]][2]
a_age = inla_spm$summary.fixed[[1]][3]
a_sex = inla_spm$summary.fixed[[1]][4]
a_bmi = inla_spm$summary.fixed[[1]][5]
b_0 = inla_spm$summary.fixed[[1]][6]
b_sex = inla_spm$summary.fixed[[1]][7]
b_bp = inla_spm$summary.fixed[[1]][8]
b_bmi = inla_spm$summary.fixed[[1]][9]
sigma_age <- bri.hyperpar.summary(inla_spm)[1]
sigma <- bri.hyperpar.summary(inla_spm)[2]
sigma2<- 0.001
c <- bri.hyperpar.summary(inla_spm)[3]
# Age effect
age_table_spm <- as.data.table(inla_spm$summary.random$m_AGE[, c("ID", "mean")])
age_table_spm[, age_round := round(ID, 2)]
age_func_spm <- gam(mean ~s(age_round),data=age_table_spm)

expl_variables_H34 <- colnames(dataH34)[which(!colnames(dataH34) %in% c("bp_4_corr", "id", "missing"))]

data_exp_H34 <- copy(dataH34[, ..expl_variables_H34])
data_exp_H34[, age_round := round(age_3, 2)]
data_exp_H34[, age_effect_spm := predict.gam(age_func_spm, data.frame(age_round = age_round)) ]


for(i in 1:n_sim){
  print(i)
  data_sim <- data_exp_H34 
  data_sim[, y_eps := rnorm(.N, 0, sigma)]
  data_sim[, y_eps_2 := rnorm(.N, 0, sigma2)]
  data_sim[, eps:= y_eps + y_eps_2]
  data_sim[, bp_4_corr := a_0 + a_bp*bp_3_corr + a_age*age_3 + a_sex*sex + a_bmi*bmi_3 + eps]
  
  data_sim[, eta_p := b_0 + b_bp*bp_3_corr + age_effect_spm + b_sex*sex + b_bp*bmi_3 + c*y_eps]
  data_sim[, p_i := exp(eta_p)/(exp(eta_p)+1)]
  data_sim[, missing := rbinom(.N, 1, p_i)]
  data_sim[, id:= 1:nH34]
  write.csv2(data_sim, file = glue::glue("/home/aurorach/data_HUNT_aurora/master/sim_data_validation/data_sim_{i}"))
  data_sim <- NULL
}


