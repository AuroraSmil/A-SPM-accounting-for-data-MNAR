#simulation_validation_mnar


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
#inla_spm <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT23.RData")
spm_sim <-readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_sim_HUNT23_naive.RData")

#number of simulations
n_sim = 100

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

# Create list of variable names
variable = c("a_0", "a_bp", "a_age",  "a_sex", "a_bmi", "b_0", "b_sex", "b_bp", "b_bmi", "sigma_age", "sigma_eps", "c")

# Store the true values of the parameters
true_values_naive_est <- data.table(variable = variable,
                                    value_true = c(inla_naive_bp$summary.fixed[[1]], 
                                                   inla_naive_m$summary.fixed[[1]], 
                                                   bri.hyperpar.summary(inla_naive_m)[1], 
                                                   bri.hyperpar.summary(inla_naive_bp)[1],
                                                   c))
age_table <- as.data.table(inla_naive_m$summary.random$m_AGE[, c("ID", "mean", "0.025quant", "0.975quant")])
age_table[, age_round := round(ID, 2)]

age_func_naive <- gam(mean ~s(age_round),data=age_table)


a_0_sim  = spm_sim$summary.fixed[[1]][1]
a_bp_sim = spm_sim$summary.fixed[[1]][2]
a_age_sim = spm_sim$summary.fixed[[1]][3]
a_sex_sim = spm_sim$summary.fixed[[1]][4]
a_bmi_sim = spm_sim$summary.fixed[[1]][5]
b_0_sim = spm_sim$summary.fixed[[1]][6]
b_sex_sim = spm_sim$summary.fixed[[1]][7]
b_bp_sim = spm_sim$summary.fixed[[1]][8]
b_bmi_sim = spm_sim$summary.fixed[[1]][9]
sigma_age_sim <- bri.hyperpar.summary(spm_sim)[1]
sigma_sim <- bri.hyperpar.summary(spm_sim)[2]
c_sim <- bri.hyperpar.summary(spm_sim)[3]
precision_age_mean_sim <- spm_sim$summary.hyperpar[[1]][1]
precision_age_sd_sim <- spm_sim$summary.hyperpar[[2]][1]
precision_eps_mean_sim <- spm_sim$summary.hyperpar[[1]][2]
precision_eps_sd_sim <- spm_sim$summary.hyperpar[[2]][2]

mean_c_sim <- spm_sim$summary.hyperpar[[1]][3]
sd_c_sim <- spm_sim$summary.hyperpar[[2]][3]

variable = c("a_0", "a_bp", "a_age",  "a_sex", "a_bmi", "b_0", "b_sex", "b_bp", "b_bmi", "sigma_age", "sigma_eps", "c", 
             "precision_age_mean", "precision_age_sd", "precision_eps_mean", "precision_eps_sd", "mean_c", "sd_c"  )
true_values_sim <- data.table(variable = variable,
                              value_true = c(spm_sim$summary.fixed[[1]], bri.hyperpar.summary(spm_sim)[1:3], 
                                             precision_age_mean_sim, precision_age_sd_sim, 
                                             precision_eps_mean_sim, precision_eps_sd_sim, mean_c_sim, sd_c_sim))




# Read data ----

dataH34<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")
#dataH34 <- dataH34[1:1000,]
n <- nrow(dataH34)
dataH34[, age_round:= round(age_3, 2)]
dataH34[, age_effect_naive := predict.gam(age_func_naive, data.frame(age_round = age_round)) ]
dataH34[, intercept := rep(b_0, n) + age_effect_naive]

expl_variables_H34 <- colnames(dataH34)[which(!colnames(dataH34) %in% c("bp_4_corr", "id", "missing"))]

data_exp_H34 <- copy(dataH34[, ..expl_variables_H34])
data_exp_H34[, age_round := round(age_3, 2)]
data_exp_H34[, age_effect_naive := predict.gam(age_func_naive, data.frame(age_round = age_round)) ]


mu_given_vec <- rep(-1, n_sim)
mu_vec <- rep(-1, n_sim)
sigma_res_given_vec <- rep(-1, n_sim)
sigma_res_vec <- rep(-1, n_sim)
# simulations starts -----
start_overall <- Sys.time()
nr_failes <- 0
i = 1
sigma2 <- 0.001
print("Starting simulatons")
while(i <=n_sim){
  print(i)
  #create data sim
  start_time <- Sys.time()
  data_sim <- data_exp_H34 # no need to copy because none of the oiginal data is overwritten
  data_sim[, y_eps := rnorm(.N, 0, sigma)]
  data_sim[, y_eps_2 := rnorm(.N, 0, sigma2)]
  data_sim[, eps:= y_eps + y_eps_2]
  data_sim[, bp_4_corr := a_0 + a_bp*bp_3_corr + a_age*age_3 + a_sex*sex + a_bmi*bmi_3 + eps]
  
  data_sim[, eta_p := b_0 + b_bp*bp_3_corr + age_effect_naive + b_sex*sex + b_bp*bmi_3 + c*y_eps]
  data_sim[, p_i := exp(eta_p)/(exp(eta_p)+1)]
  data_sim[, missing := rbinom(.N, 1, p_i)]
  data_sim[, id:= 1:n]
  
  joint_inla <- pred_spm_given_m(data_sim, true_values_sim,n)
  #summary(joint_inla)
  joint_inla_n <- pred_spm_bp(data_sim, true_values_sim, n)
  #summary(joint_inla_n)
  
  res_given <- data_sim$bp_4_corr- joint_inla$summary.fitted.values$mean[1:n]
  res_n <- data_sim$bp_4_corr- joint_inla_n$summary.fitted.values$mean[1:n]
  
  mu_given<- mean(abs(res_given))
  mu<- mean(abs(res_n))
  print(mu_given)
  print(mu)
  
  data_sim <- NULL
  #sigma_given <- sqrt(sum(na.omit(res_given - mu_given)^2)/n^2)
  #sigma <- sqrt(sum(na.omit(res_n - mu)^2)/n^2)
  
  mu_given_vec[i] <- mu_given
  mu_vec[i] <- mu
  #sigma_res_given_vec[i] <- sigma_given
  #sigma_res_vec[i] <- sigma
  
  i <- i + 1
  
}

today <- today()
if(!dir.exists(glue("results/simulations/{today}"))){
  dir.create(glue("results/simulations/{today}"))
}
if(!dir.exists(glue("results/simulations/{today}/1"))){
  dir.create(glue("results/simulations/{today}/1"))
}

dir <- glue("results/simulations/{today}/1")
write.csv(mu_given_vec, file = glue("{dir}/{n_sim}_naive_mu_given_vec"), row.names = FALSE)
write.csv(mu_vec, file = glue("{dir}/{n_sim}_naive_mu_vec"), row.names = FALSE)
#write.csv(sigma_res_given_vec, file = glue("{dir}/{n_sim}_sigma_res_given_vec"), row.names = FALSE)
#write.csv(sigma_res_vec, file = glue("{dir}/{n_sim}_sigma_res_vec"), row.names = FALSE)
#mu_vec_given <- rbind(mu_given_vec)
#mu_vec <- rbind(mu_vec)
# res_data <- data.table(mu_given_vec = as.list(mu_given_vec), 
#                        mu_vec = mu_vec, 
#                        sigma_res_given_vec = sigma_res_given_vec, 
#                        sigma_res_vec = sigma_res_vec)
# 
# write.csv2(res_data, file = glue("{dir}/{n_sim}_res_data"), row.names = FALSE)