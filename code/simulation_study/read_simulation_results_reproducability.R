library(data.table)
library(ggplot2)

library(cowplot)
library(gridExtra)
library(grid)
library(gridtext)
library(latex2exp)

source("code/simulation_study/model_functions.R")

dir_spm_mnar <- "results/simulations/2021-12-17"
dir_spm_mar <- "results/simulations/2021-12-18"
dir_naive_mnar <- "results/simulations/2021-12-18"
dir_naive_mar <- "results/simulations/2021-12-19"
### SPM estimates MNAR----

spm_est_coverage_spm_mnar<-  as.data.table(read.csv(file = glue::glue("{dir_spm_mnar}/SPM_parameters_MNAR/100_coverage_spm")))
spm_est_coverage_naive_mnar<-  as.data.table(read.csv(file = glue::glue("{dir_spm_mnar}/SPM_parameters_MNAR/100_coverage_naive")))

spm_est_parameter_est_spm_mnar <- as.data.table(read.csv(file = glue::glue("{dir_spm_mnar}/SPM_parameters_MNAR/100_parameter_est_spm")))
spm_est_parameter_est_naive_mnar <-  as.data.table(read.csv(file = glue::glue("{dir_spm_mnar}/SPM_parameters_MNAR/100_parameter_est_naive")))

spm_est_true_values_mnar <- as.data.table(read.csv(file = glue::glue("{dir_spm_mnar}/SPM_parameters_MNAR/true_values")))

spm_est_mnar <- create_sim_data(spm_est_parameter_est_spm_mnar, 
                          spm_est_coverage_spm_mnar, 
                          spm_est_parameter_est_naive_mnar, 
                          spm_est_coverage_naive_mnar,
                          spm_est_true_values_mnar)

spm_est_mnar$simulation_results_sumary

q_spm_est_mnar<- plot_sim_results(spm_est_mnar$simulation_results, spm_est_mnar$simulation_results_sumary )
ggsave("images/simulation_study_results_MNAR_spm_est_100.pdf",
       q_spm_est_mnar,
       width = 15,
       height = 20,
       units = "cm")

setnames(spm_est_mnar$simulation_results_sumary, 
         c("true_value", "mean_spm", "bias_spm", "coverage_spm", 
           "mean_naive", "bias_naive", "coverage_naive", "diff_bias"),
         c("True value", "Mean", "Bias", "Coverage", "Mean", "Bias", "Coverage", "Difference in bias"))
spm_est_mnar$simulation_results_sumary[, variable:= c("$alpha_0$", "$alpha_{BP}$", "$alpha_{age}$", "$alpha_{sex}$", "$alpha_{bmi}$",
                                                      "$beta_0$", "$beta_{sex}$", "$beta_{BP}$", "$beta_{BMI}$", "$sigma_{age}$", "$sigma_{epsilon}$", "c")]
spm_est_mnar$simulation_results_sumary[, row_order:= c(1,2,3,5,4,6,8,9,7,10,11,12)]
spm_est_mnar$simulation_results_sumary<- spm_est_mnar$simulation_results_sumary[order(row_order)]
spm_est_mnar$simulation_results_sumary[,row_order:= NULL]
length(which(is.na(spm_est_parameter_est_spm_mnar$sigma_eps)))
format(spm_est_mnar$simulation_results_sumary, digits = 2)
write.csv(format(spm_est_mnar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_spm_param_est_mnar", row.names = FALSE)

### MAR spm estimates ----

spm_est_coverage_spm_mar<-  as.data.table(read.csv(file = glue::glue("{dir_spm_mar}/SPM_parameters_MAR/100_coverage_spm")))
spm_est_coverage_naive_mar<-  as.data.table(read.csv(file = glue::glue("{dir_spm_mar}/SPM_parameters_MAR/100_coverage_naive")))

spm_est_parameter_est_spm_mar <- as.data.table(read.csv(file = glue::glue("{dir_spm_mar}/SPM_parameters_MAR/100_parameter_est_spm")))
spm_est_parameter_est_naive_mar <-  as.data.table(read.csv(file = glue::glue("{dir_spm_mar}/SPM_parameters_MAR/100_parameter_est_naive")))


spm_est_true_values_mar <- as.data.table(read.csv(file = glue::glue("{dir_spm_mar}/SPM_parameters_MAR/true_values")))
spm_est_mar <- create_sim_data(spm_est_parameter_est_spm_mar, 
                                spm_est_coverage_spm_mar, 
                                spm_est_parameter_est_naive_mar, 
                                spm_est_coverage_naive_mar,
                                spm_est_true_values_mar)

spm_est_mar$simulation_results_sumary

q_spm_est_mar<- plot_sim_results(spm_est_mar$simulation_results, spm_est_mar$simulation_results_sumary )

ggsave("images/simulation_study_results_MAR_spm_est_100.pdf",
       q_spm_est_mar,
       width = 15,
       height = 20,
       units = "cm")

setnames(spm_est_mar$simulation_results_sumary, 
         c("true_value", "mean_spm", "bias_spm", "coverage_spm", 
           "mean_naive", "bias_naive", "coverage_naive", "diff_bias"),
         c("True value", "Mean", "Bias", "Coverage", "Mean", "Bias", "Coverage", "Difference in bias"))
spm_est_mar$simulation_results_sumary[, variable:= c("$alpha_0$", "$alpha_{BP}$", "$alpha_{age}$", "$alpha_{sex}$", "$alpha_{bmi}$",
                                                     "$beta_0$", "$beta_{sex}$", "$beta_{BP}$", "$beta_{BMI}$", "$sigma_{age}$", "$sigma_{epsilon}$", "c")]
spm_est_mar$simulation_results_sumary[, row_order:= c(1,2,3,5,4,6,8,9,7,10,11,12)]
spm_est_mar$simulation_results_sumary<- spm_est_mar$simulation_results_sumary[order(row_order)]
spm_est_mar$simulation_results_sumary[,row_order:= NULL]
length(which(is.na(spm_est_parameter_est_spm_mar$sigma_eps)))

format(spm_est_mar$simulation_results_sumary, digits = 2)
write.csv(format(spm_est_mar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_spm_param_est_mar", row.names = FALSE)

######################################

##################################

### MNAR naive estimates ----
naive_est_coverage_spm_mnar<-  as.data.table(read.csv(file = glue::glue("{dir_naive_mnar}/naive_parameters_MNAR/100_coverage_spm")))
naive_est_coverage_naive_mnar<-  as.data.table(read.csv(file = glue::glue("{dir_naive_mnar}/naive_parameters_MNAR/100_coverage_naive")))

naive_est_parameter_est_spm_mnar <- as.data.table(read.csv(file = glue::glue("{dir_naive_mnar}/naive_parameters_MNAR/100_parameter_est_spm")))
naive_est_parameter_est_naive_mnar <-  as.data.table(read.csv(file = glue::glue("{dir_naive_mnar}/naive_parameters_MNAR/100_parameter_est_naive")))

naive_est_true_values_mnar <- as.data.table(read.csv(file = glue::glue("{dir_naive_mnar}/naive_parameters_MNAR/true_values")))


naive_est_mnar <- create_sim_data(naive_est_parameter_est_spm_mnar, 
                               naive_est_coverage_spm_mnar, 
                               naive_est_parameter_est_naive_mnar, 
                               naive_est_coverage_naive_mnar,
                               naive_est_true_values_mnar)

naive_est_mnar$simulation_results_sumary

q_naive_est_mnar<- plot_sim_results(naive_est_mnar$simulation_results, naive_est_mnar$simulation_results_sumary )

ggsave("images/simulation_study_results_MNAR_naive_est_100.pdf",
       q_naive_est_mnar,
       width = 15,
       height = 20,
       units = "cm")

setnames(naive_est_mnar$simulation_results_sumary, 
         c("true_value", "mean_spm", "bias_spm", "coverage_spm", 
           "mean_naive", "bias_naive", "coverage_naive", "diff_bias"),
         c("True value", "Mean", "Bias", "Coverage", "Mean", "Bias", "Coverage", "Difference in bias"))
naive_est_mnar$simulation_results_sumary[, variable:= c("$alpha_0$", "$alpha_{BP}$", "$alpha_{age}$", "$alpha_{sex}$", "$alpha_{bmi}$",
                                                        "$beta_0$", "$beta_{sex}$", "$beta_{BP}$", "$beta_{BMI}$", "$sigma_{age}$", "$sigma_{epsilon}$", "c")]
naive_est_mnar$simulation_results_sumary[, row_order:= c(1,2,3,5,4,6,8,9,7,10,11,12)]
naive_est_mnar$simulation_results_sumary <- naive_est_mnar$simulation_results_sumary[order(row_order)]
naive_est_mnar$simulation_results_sumary[,row_order:= NULL]
length(which(is.na(naive_est_parameter_est_spm_mnar$sigma_eps)))

format(naive_est_mnar$simulation_results_sumary, digits = 2)
write.csv(format(naive_est_mnar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_naive_est_mnar", row.names = FALSE)

### MAR naive estimates ----

naive_est_coverage_spm_mar<-  as.data.table(read.csv(file = glue::glue("{dir_naive_mar}/naive_parameters_MAR/100_coverage_spm")))
naive_est_coverage_naive_mar<-  as.data.table(read.csv(file = glue::glue("{dir_naive_mar}/naive_parameters_MAR/100_coverage_naive")))

naive_est_parameter_est_spm_mar <- as.data.table(read.csv(file = glue::glue("{dir_naive_mar}/naive_parameters_MAR/100_parameter_est_spm")))
naive_est_parameter_est_naive_mar <-  as.data.table(read.csv(file = glue::glue("{dir_naive_mar}/naive_parameters_MAR/100_parameter_est_naive")))

naive_est_true_values_mar <- as.data.table(read.csv(file = glue::glue("{dir_naive_mar}/naive_parameters_MAR/true_values")))


naive_est_mar <- create_sim_data(naive_est_parameter_est_spm_mar, 
                                  naive_est_coverage_spm_mar, 
                                  naive_est_parameter_est_naive_mar, 
                                  naive_est_coverage_naive_mar,
                                  naive_est_true_values_mar)

naive_est_mar$simulation_results_sumary

q_naive_est_mar<- plot_sim_results(naive_est_mar$simulation_results, naive_est_mar$simulation_results_sumary )


ggsave("images/simulation_study_results_MAR_naive_est_100.pdf",
       q_naive_est_mar,
       width = 15,
       height = 20,
       units = "cm")

setnames(naive_est_mar$simulation_results_sumary, 
         c("true_value", "mean_spm", "bias_spm", "coverage_spm", 
           "mean_naive", "bias_naive", "coverage_naive", "diff_bias"),
         c("True value", "Mean", "Bias", "Coverage", "Mean", "Bias", "Coverage", "Difference in bias"))
naive_est_mar$simulation_results_sumary[, variable:= c("$alpha_0$", "$alpha_{BP}$", "$alpha_{age}$", "$alpha_{sex}$", "$alpha_{bmi}$",
                                                        "$beta_0$", "$beta_{sex}$", "$beta_{BP}$", "$beta_{BMI}$", "$sigma_{age}$", "$sigma_{epsilon}$", "c")]
naive_est_mar$simulation_results_sumary[, row_order:= c(1,2,3,5,4,6,8,9,7,10,11,12)]
naive_est_mar$simulation_results_sumary<- naive_est_mar$simulation_results_sumary[order(row_order)]
naive_est_mar$simulation_results_sumary[,row_order:= NULL]
length(which(is.na(naive_est_parameter_est_spm_mar$sigma_eps)))

format(naive_est_mar$simulation_results_sumary, digits = 2)
write.csv(format(naive_est_mar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_naive_est_mar", row.names = FALSE)





