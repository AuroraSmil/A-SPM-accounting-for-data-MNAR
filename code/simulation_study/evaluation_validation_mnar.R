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


dir <- "results/simulations/2022-01-05/1/100"
mu_given_vec<- read.csv(file = glue("{dir}_naive_mu_given_vec"))
mu_vec<- read.csv(file = glue("{dir}_naive_mu_vec"))
#write.csv(sigma_res_given_vec, file = glue("{dir}/{n_sim}_sigma_res_given_vec"), row.names = FALSE)
#write.csv(sigma_res_vec, file = glue("{dir}/{n_sim}_sigma_res_vec"), row.names = FALSE)
#mu_vec_given <- rbind(mu_given_vec)
#mu_vec <- rbind(mu_vec)
res_data <- data.table(mu_given_vec = mu_given_vec[1:100,1], 
                       mu_vec = mu_vec[1:100,1])
res_data[, diff := mu_given_vec- mu_vec]
q <- ggplot(data = res_data )
q <- q + geom_density(aes(x = mu_given_vec, colour = "m known", linetype = "m known"))
q <- q + geom_density(aes(x = mu_vec, colour = "m unknown", linetype = "m unknown"))
q <- q + xlab("MAE") #+ xlim(0.046,0.0465)
q <- q + scale_color_discrete(labels = expression(hat(BP[F])~"|"~m,hat(BP[F])))
q <- q + scale_linetype_discrete(labels = expression(hat(BP[F])~"|"~m,hat(BP[F])))
q <- q + labs(color  = "Guide name", linetype = "Guide name") + ylab("Density")
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q
ggsave("images/validation_mnar_naive_param.pdf",
       q,
       width = 17,
       height = 10,
       units = "cm")

# q <- ggplot(data = res_data )
# q <- q + geom_density(aes(x = diff))
# q <- q + xlab(TeX("(\\bar{|True $BP_F$ - predicted $\\bar{BP_F|m}$|})-\\bar{(|True $BP_F$ - predicted $\\bar{BP_F}$|})")) + xlim(-0.0017,0.0005)
# q <- q + geom_vline(xintercept = 0, linetype = "dotted")
# q <- q+ ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
# q
# 
# 
# ggsave("images/validation_mnar_diff_naive_param.pdf",
#        q,
#        width = 17,
#        height = 10,
#        units = "cm")

dir <- "results/simulations/2022-01-05/100"
mu_given_vec_spm<- read.csv(file = glue("{dir}_mu_given_vec"))
mu_vec_spm<- read.csv(file = glue("{dir}_mu_vec"))
res_data_spm <- data.table(mu_given_vec = mu_given_vec_spm[1:100,1], 
                           mu_vec = mu_vec_spm[1:100,1])
res_data_spm[, diff := mu_given_vec- mu_vec]
q <- ggplot(data = res_data_spm )
q <- q + geom_density(aes(x = mu_given_vec, colour = "m known", linetype = "m known"))
q <- q + geom_density(aes(x = mu_vec, colour = "m unknown", linetype = "m unknown"))
q <- q + xlab("MAE") #+ xlim(0.046,0.0465)
q <- q + scale_color_discrete(labels = expression(hat(BP[F])~"|"~m,hat(BP[F])))
q <- q + scale_linetype_discrete(labels = expression(hat(BP[F])~"|"~m,hat(BP[F])))
q <- q + labs(color  = "Guide name", linetype = "Guide name") + ylab("Density") 
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q


ggsave("images/validation_mnar_spm_param.pdf",
       q,
       width = 17,
       height = 10,
       units = "cm")

# q <- ggplot(data = res_data_spm )
# q <- q + geom_density(aes(x = diff))
# q <- q + xlab(TeX("(\\bar{|True $BP_F$ - predicted $\\bar{BP_F|m}$|})-\\bar{(|True $BP_F$ - predicted $\\bar{BP_F}$|})")) + xlim(-0.0017,0.0005)
# q <- q + geom_vline(xintercept = 0, linetype = "dotted")
# q <- q+ ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
# q
# 
# 
# ggsave("images/validation_mnar_diff_spm_param.pdf",
#        q,
#        width = 17,
#        height = 10,
#        units = "cm")


q <- ggplot(data = res_data_spm )
q <- q + geom_density(aes(x = diff, colour = "Data MNAR", linetype = "Data MNAR"))
q <- q + geom_density(data= res_data, aes(x = diff, colour = "Data MAR", linetype = "Data MAR"))
q <- q + xlab(TeX("$MAE(\\hat{BP_F}|m)- MAE(\\hat{BP_F})$")) + xlim(-0.0017,0.0005)

#q <- q + xlab(TeX("(\\bar{|True $BP_F$ - predicted $\\bar{BP_F|m}$|})-\\bar{(|True $BP_F$ - predicted $\\bar{BP_F}$|})")) + xlim(-0.0017,0.0005)
q <- q + geom_vline(xintercept = 0, linetype = "dotted")
q <- q+ ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + labs(colour  = "Guide name", linetype = "Guide name")
q

ggsave("images/validation_mnar_diff_param.pdf",
       q,
       width = 17,
       height = 10,
       units = "cm")

