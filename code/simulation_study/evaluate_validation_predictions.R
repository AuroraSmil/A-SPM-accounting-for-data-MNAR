# Eval validataion predictions. 


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
today <- today()
dir <- glue("results/simulations/2022-01-18")
n_sim <- 100
crps_n <- as.data.table(read.csv(file = glue("{dir}/naive/{n_sim}_crps_sim_model")))
crps_spm <-  as.data.table(read.csv(file = glue("{dir}/SPM/{n_sim}_crps_sim_model")))
brier_n <-  as.data.table(read.csv(file = glue("{dir}/naive/{n_sim}_brier_sim_model")))
brier_spm <-  as.data.table(read.csv(file = glue("{dir}/SPM/{n_sim}_brier_sim_model")))

setnames(crps_n, "crps", "crps_n")
setnames(crps_spm, "crps", "crps_spm")
setnames(brier_n, "brier", "brier_n")
setnames(brier_spm, "brier", "brier_spm")

crps <- merge(crps_n, crps_spm, on= c(missing, sim_id))
brier<- merge(brier_n, brier_spm, on = c(missing, sim_id))

crps[, crps_spm_minus_n := crps_spm-crps_n]

crps_agg <- crps[, .(crps_spm_mean = mean(crps_spm),
                           crps_n_mean = mean(crps_n),
                           crps_spm_minus_n_mean = mean(crps_spm_minus_n)),
                       keyby = .(missing)]

brier[, brier_spm_minus_n := brier_spm - brier_n]
brier_agg <- brier[, .(brier_spm_mean = mean(brier_spm),
                            brier_n_mean = mean(brier_n),
                            brier_spm_minus_n_mean = mean(brier_spm_minus_n)),
                       keyby = .(missing)]

# Plot the final results
q <- ggplot(data = crps, aes(group = as.factor(missing), colour = as.factor(missing), linetype = as.factor(missing)))
q <- q + geom_density(aes(x = crps_spm_minus_n))
q <- q + xlab(TeX("$\\bar{CRPS}_{SPM}- \\bar{CRPS}_{naive}$"))
q <- q + scale_color_discrete(labels = c("All", "Pressent", "Missing"))
q <- q + scale_linetype_discrete(labels = c("All", "Pressent", "Missing"))
q <- q + labs(color  = "Guide name", linetype = "Guide name") + ylab("Density") 
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q
ggsave("images/diff_crps.pdf",
       q,
       width = 17,
       height = 10,
       units = "cm")

q <- ggplot(data = brier, aes(group = as.factor(missing), colour = as.factor(missing), linetype = as.factor(missing)))
q <- q + geom_density(aes(x = brier_spm_minus_n))
q <- q + xlab(TeX("$Brier_{SPM}- \\Brier_{naive}$"))
q <- q+ ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + scale_color_discrete(labels = c("All", "Pressent", "Missing"))
q <- q + scale_linetype_discrete(labels = c("All", "Pressent", "Missing"))
q <- q + labs(color  = "Guide name", linetype = "Guide name") + ylab("Density") 
q
ggsave("images/diff_brier.pdf",
       q,
       width = 17,
       height = 10,
       units = "cm")

write.csv(crps_agg, file = glue("{dir}/{n_sim}_crps_agg_param_est_spm_mnar"), row.names = FALSE)
write.csv(brier_agg, file = glue("{dir}/{n_sim}_brier_agg_param_est_spm_mnar"), row.names = FALSE)

