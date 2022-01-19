library(INLA)
library(brinla)
#INLA:::inla.dynload.workaround()
library(tidyverse)
library(data.table)
library(latex2exp)
library(tikzDevice)
library(mgcv)
source("code/validation/validation_func.R")

# read SPM and model parameters ####
inla_spm <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT23.RData")

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
c <- bri.hyperpar.summary(inla_spm)[3]

precision_age_mean <- inla_spm$summary.hyperpar[[1]][1]
precision_age_sd <- inla_spm$summary.hyperpar[[2]][1]
precision_eps_mean <- inla_spm$summary.hyperpar[[1]][2]
precision_eps_sd <- inla_spm$summary.hyperpar[[2]][2]

mean_c <- inla_spm$summary.hyperpar[[1]][3]
sd_c <- inla_spm$summary.hyperpar[[2]][3]


variable = c("a_0", "a_bp", "a_age",  "a_sex", "a_bmi", "b_0", "b_sex", "b_bp", "b_bmi", "sigma_age", "sigma_eps", "c", 
             "precision_age_mean", "precision_age_sd", "precision_eps_mean", "precision_eps_sd", "mean_c", "sd_c"  )
true_values <- data.table(variable = variable,
                          value_true = c(inla_spm$summary.fixed[[1]], bri.hyperpar.summary(inla_spm)[1:3], 
                                         precision_age_mean, precision_age_sd, precision_eps_mean, precision_eps_sd, mean_c, sd_c))

# Read data ---- 
data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")
data_small <- as.data.table(data)#[1:1000])
n <- nrow(data_small)
data_sim <- data_small

age_table <- as.data.table(inla_spm$summary.random$m_AGE[, c("ID", "mean", "0.025quant", "0.975quant")])
age_table[, age_round := round(ID, 2)]
data_sim[, age_round:= round(age_3, 2)]

age_func_spm <- gam(mean ~s(age_round),data=age_table)
data_sim[, age_effect_spm := predict.gam(age_func_spm, data.frame(age_round = age_round)) ]
data_sim[, intercept := rep(b_0, n) + age_effect_spm]

joint_inla <- pred_spm_given_m(data_sim, true_values,n)

summary(joint_inla)
length(joint_inla$summary.fitted.values$mean)

joint_inla_n <- pred_spm_bp(sim_data, true_values, n)


# Comparing residuals -----
res_given <- data_sim$bp_4_corr- joint_inla$summary.fitted.values$mean[1:n]
res_n <- data_sim$bp_4_corr- joint_inla_n$summary.fitted.values$mean[1:n]

diff_abs_res_given_minus_res_n <- abs(res_given) - abs(res_n)
mean(na.omit(diff_abs_res_given_minus_res_n))
diff_res_given_minus_res_n <- res_given - res_n
mean(abs(na.omit(diff_res_given_minus_res_n)))
mean((na.omit(diff_res_given_minus_res_n)))

mu_given<- mean(abs(na.omit(res_given)))
mu<- mean(abs(na.omit(res_n)))

mean(abs(na.omit(res_given^2)))
mean(abs(na.omit(res_n^2)))

sigma_given <- sqrt(sum(na.omit(res_given - mu_given)^2)/n^2)
sigma <- sqrt(sum(na.omit(res_n - mu)^2)/n^2)
mu_given
mu
mu_given - mu
sigma_given
sigma
mu_given - sigma_given
mu + sigma
# ##### Small hypothesis in significanse
# 
# # H0: the two models preform equally H1 missing status enhances the predictions
# #Under H0 al residuals come from the same distribution hence we can mix them up
# 
# significantse_cut <- mean(abs(na.omit(res_n))) -  mean(abs(na.omit(res_given)))
# all_res <- na.omit(c(res_given, res_n))
# count <- 0
# n_hyp <- 10000
# for(i in 1:n_hyp){
#   n <- length(all_res)
#   index <- sample(1:n, n/2)
#   temp_given <- all_res[index]
#   
#   temp <- all_res[-index]
#   diff <- mean(abs(na.omit(temp))) -  mean(abs(na.omit(temp_given)))
#   if (diff >significantse_cut){
#     count <- count + 1
#   }
# }
# 
# count/n_hyp
# 
# 
# # Plot resiguals
# data_res <- data.table(res_given = res_given, 
#                        res_n = res_n,
#                        bp_4_corr = data_sim$bp_4_corr, 
#                        bp_4_spm = joint_inla$summary.fitted.values$mean[1:n],
#                        bp_4_n = joint_inla_n$summary.fitted.values$mean[1:n], 
#                        missing = as_factor(data_sim$missing))
# q <- ggplot(data= data_res, aes(x =bp_4_spm, y = res_given, colour = "spm"))
# q <- q + geom_point() 
# q <- q + geom_point( aes(x =bp_4_n, y = res_n, colour = "naive"))
# q <- q + geom_smooth(method = "lm", se= FALSE)
# q <- q + geom_smooth(aes(x =bp_4_n, y = res_n, colour = "naive"), method = "lm", se = FALSE)
# q <- q + geom_hline(yintercept = 0, color = "red")
# q
# 
# 
# q <- ggplot(data_res, aes(x = bp_4_spm, group = missing, colour = "information about missing status", linetype = missing))
# q <- q + geom_density()
# q <- q + geom_density(aes(x = bp_4_n, group = missing, colour = "no information", linetype = missing))
# q <- q + geom_density(aes(x = bp_4_corr, group = missing, colour = "true values", linetype = missing))
# q <- q + geom_vline(xintercept = mean(na.omit(data_res$bp_4_corr)))
# q <- q + geom_vline(xintercept = mean(na.omit(data_res$bp_4_spm)), colour = "red")
# q <- q + geom_vline(xintercept = mean(na.omit(data_res$bp_4_n)), colour = "blue")
# q



