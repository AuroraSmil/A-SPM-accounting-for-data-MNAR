library(INLA)
library(brinla)
#INLA:::inla.dynload.workaround()
library(tidyverse)
library(data.table)
library(latex2exp)
library(tikzDevice)
library(lubridate)
library(glue)
library(lubridate)

source("code/simulation_study/model_functions.R")

# Read data
data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")#[1:1000]

rel_col_fixed <- c("mean", "0.025quant", "0.975quant")
rel_col_hyper <- c("mean", "q0.025", "q0.975")
parameters = c("a_0", "a_bp", "a_age",  "a_sex", "a_bmi", "b_0", "b_sex", "b_bp", "b_bmi", "sigma_age", "sigma_eps", "c")

# plot_c <- function(spm_model){
#   c<- spm_model$marginals.hyperpar$`Beta for m_eps1`
#   
#   q <- ggplot(data = as.data.frame(c), aes(x = x))
#   q <- q + geom_line(aes(y = y, colour = "SPM"), linetype = "dashed")
#   q <- q + xlab(unname(TeX(c("c"))))
#   q <- q + ylab("Density")
#   q <- q + ylab(element_blank())
#   q <- q + theme(legend.title = element_blank())
#   q <- q + scale_color_manual(values=c("#00BFC4"))
#   return(q)
# }
# 
# plot_age_effect <- function(spm_model){
#   data_random_age_inla <-  spm_model$summary.random$m_AGE
#   q <- ggplot(data = data_random_age_inla, aes(x = ID))
#   q <- q + geom_line(aes(y = mean, colour = "SPM", linetype = "SPM"))
#   q <- q + geom_ribbon(aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "SPM"), alpha = 0.5)
#   q <- q + labs(colour  = "Guide name", linetype = "Guide name")
#   q <- q + scale_linetype_manual(values = c("solid"))
#   q_age <- q + xlab(TeX("Age")) + ylab("Non-linear effect") +guides(fill = FALSE) + theme( legend.position = "bottom", legend.title = element_blank())
#   return(q_age)
# }

# Priors for c ---- 
priors <- data.table(
  sigma = c(1, 1,10,100, 1, 10, 100, 100),
  mean = c(10, 0,0,0,1,1,1,10)
)

# priors <- data.table(
#   sigma = c(1), 
#   mean = c(100)
# )
priors[, log_prec:= log(1/sigma^2)]
priors[, prec := 1/sigma^2]
priors[, model:= paste0("mu = ", mean," sigma = ",sigma)]

################# Run simulations ################

#q_age <- vector("list", length = nrow(priors))
for( i in 1:nrow(priors)){
  print(glue::glue("Round: {i}. mean: {priors[i]$mean}. precision: {priors[i]$prec}."))
  spm <- run_spm(data, priors[i]$mean, priors[i]$prec, verbose = F)
  
  #q_age[[i]]<- plot_age_effect(spm)
  
  print(spm$summary.fixed)
  print(bri.hyperpar.summary(spm))
  
  
  summary_fixed <- as.data.table(spm$summary.fixed)[, ..rel_col_fixed]
  summary_hyper <- as.data.table(bri.hyperpar.summary(spm))[, ..rel_col_hyper]
  
  summary <- rbindlist(list(summary_fixed, summary_hyper ), use.names = FALSE)
  
  summary[, parameter := parameters]
  summary[, model := priors[i]$model]
  summary[, model_mean := priors[i]$mean]
  summary[, sigma := priors[i]$sigma]
  summary_age_effekt_mean <- as.data.table(spm$summary.random$m_AGE)[, c("ID", "mean")]
  summary_age_effekt_mean[, model := priors[i]$model]
  summary_age_effekt_mean[, model_mean := priors[i]$mean]
  summary_age_effekt_mean[, sigma := priors[i]$sigma]
  
  if (i == 1){
    summary_total<- summary
    summary_age_effekt_mean_total <- summary_age_effekt_mean
  }else{
    summary_total <- rbindlist(list(summary_total, summary ))
    summary_age_effekt_mean_total <-rbindlist(list(summary_age_effekt_mean_total, summary_age_effekt_mean))
  }
}

today <- today()
write.csv(summary_total, file = glue::glue("results/prior_sensitivity/priors_{today}"), row.names = FALSE)
write.csv(summary_age_effekt_mean_total, file = glue::glue("results/prior_sensitivity/age_effect_{today}"), row.names = FALSE)

# Plots 
variable_names <- c(
  `a_0` = "alpha[0]",
  `a_bp` = "alpha[BP]",
  `a_age` = "alpha[age]",
  `a_sex` = "alpha[sex]",
  `a_bmi` = "alpha[BMI]", 
  `b_0` = "beta[0]",
  `b_sex` = "beta[sex]",
  `b_bp` = "beta[BP]",
  `b_bmi` = "beta[BMI]",
  `sigma_age` = "sigma[age]",
  `sigma_eps` = "sigma[epsilon]",
  `c`  = "c"    
)


summary_total[, model_name := TeX(paste0("\\mu = ",model_mean, " \\sigma = ", sigma))]
q <- ggplot(data = summary_total, aes( x = mean, y = model))
q <- q + geom_point(aes(x = mean))
q <- q + geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`), height= .2)
q <- q + ylab("Models") + xlab(element_blank())
q <- q + scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
q<- q + facet_wrap(~parameter, scales = "free_x", labeller = labeller(parameter = as_labeller(variable_names, label_parsed)))
q

ggsave("images/prior_sensitivity_c.pdf",
       q,
       width = 20,
       height = 15,
       units = "cm")

q <- ggplot(data = summary_age_effekt_mean_total, aes( x = ID, y = mean, group = model, colour = model, linetype = model))
q <- q + geom_line()
q <- q + ylab("Mean age effect") + xlab("Age")
q <- q+ theme( legend.position = "bottom", legend.title = element_blank())
q

ggsave(glue::glue("images/prior_sensitivity_c_age_effect.pdf"),
       q,
       width = 20,
       height = 10,
       units = "cm")

