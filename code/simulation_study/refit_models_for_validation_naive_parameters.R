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

# Read models and data ---- 


dataH23 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")

nH23 = nrow(dataH23)
# Store spm parameters 

inla_naive_bp<-readRDS(file = "/home/aurorach/data_HUNT_aurora/master/naive_bp_inla_HUNT2.RData")
inla_naive_m<-readRDS(file = "/home/aurorach/data_HUNT_aurora/master/naive_m_inla_HUNT2.RData")

a_0  = inla_naive_bp$summary.fixed[[1]][1]
a_bp = inla_naive_bp$summary.fixed[[1]][2]
a_age = inla_naive_bp$summary.fixed[[1]][3]
a_sex = inla_naive_bp$summary.fixed[[1]][4]
a_bmi = inla_naive_bp$summary.fixed[[1]][5]
b_0 = inla_naive_m$summary.fixed[[1]][1]
b_sex = inla_naive_m$summary.fixed[[1]][2]
b_bp = inla_naive_m$summary.fixed[[1]][3]
b_bmi = inla_naive_m$summary.fixed[[1]][4]


sigma_age <- bri.hyperpar.summary(inla_naive_m)[1]
sigma <- bri.hyperpar.summary(inla_naive_bp)[1]
c <- 0

sigma2 <- 0.001
# Age effect
age_table_naive <- as.data.table(inla_naive_m$summary.random$m_AGE[, c("ID", "mean")])
age_table_naive[, age_round := round(ID, 2)]
age_func_naive <- gam(mean ~s(age_round),data=age_table_naive)


# Explanaotry variables 
expl_variables_H23 <- colnames(dataH23)[which(!colnames(dataH23) %in% c("bp_3_corr", "id", "missing"))]

data_exp_H23 <- copy(dataH23[, ..expl_variables_H23])
data_exp_H23[, age_round := round(age_2, 2)]
data_exp_H23[, age_effect_naive := predict.gam(age_func_naive, data.frame(age_round = age_round)) ]


# Simulated data to refit SPM and the naive model ---- 

data_sim <- data_exp_H23 # no need to copy because none of the oiginal data is overwritten
data_sim[, y_eps := rnorm(.N, 0, sigma)]
data_sim[, y_eps_2 := rnorm(.N, 0, sigma2)]
data_sim[, eps:= y_eps + y_eps_2]
data_sim[, bp_3_corr := a_0 + a_bp*bp_2_corr + a_age*age_2 + a_sex*sex + a_bmi*bmi_2 + eps]

data_sim[, eta_p := b_0 + b_bp*bp_2_corr + age_effect_naive + b_sex*sex + b_bp*bmi_2 + c*y_eps]
data_sim[, p_i := exp(eta_p)/(exp(eta_p)+1)]
data_sim[, missing := rbinom(.N, 1, p_i)]
data_sim[missing == 1, bp_3_corr := NA]
data_sim[, id:= 1:nH23]

saveRDS(data_sim, file = "/home/aurorach/data_HUNT_aurora/master/HUNT23_sim_data_naive.Rdata")


############# start here ###############3

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



data_sim <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/HUNT23_sim_data.Rdata")

# naive_m_sim <- run_naive_m(data_sim, compute_all = TRUE)
# saveRDS(naive_m_sim, file = "/home/aurorach/data_HUNT_aurora/master/inla_naive_m_sim_HUNT23.RData")
# naive_bp_sim <- run_naive_bp(data_sim, compute_all = TRUE)
# saveRDS(naive_bp_sim, file = "/home/aurorach/data_HUNT_aurora/master/inla_naive_bp_sim_HUNT23.RData")
spm_sim <- run_spm(data_sim, verbose = F, compute_all = TRUE)
saveRDS(spm_sim, file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_sim_HUNT23_naive.RData")


