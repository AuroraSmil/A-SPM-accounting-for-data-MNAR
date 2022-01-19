# 
# library(INLA)
# #library(inlabru)
# library(data.table)
# library(brinla)
# library(ggplot2)
# library(dplyr)
# library(DescTools)
# library(latex2exp)
# library(tikzDevice)
# library(scoringRules)
# library(glue)
# library(mgcv)


assosiatoin_effect <- function(dir){
  # Read model and data ----
  inla_spm <- readRDS(file = glue::glue("{dir}spm_inla_HUNT23.RData"))
  inla_spm$summary.fixed[6:9,]
  
  # Get age effect
  age_table <- as.data.table(inla_spm$summary.random$m_AGE[, c("ID", "mean", "0.025quant", "0.975quant")])
  age_table[, age_round := round(ID, 2)]
  age_func <- gam(mean ~s(age_round),data=age_table)
  
  ##### Simulated individuals ----
  n_samples <- 1000
  set.seed(123)
  data_fake <- data.table(
    id = c(1,2,3),
    bp_2_corr = c(-1.7, 0, 1.7),
    age_round = c(-0.50, 0.00, 0.50),
    age_2 = c(-0.5, 0, 0.5),
    bmi_2 = c(-0.8, 0, 0.8),
    sex = c(0,0,0),
    name = c("Current systolic bloodpressure = -1.7 ", "Current systolic bloodpressure = 0",
             "Current systolic bloodpressure = 1.7")
  )
  
  data_fake <- data.table(
    id = c(1,2,3),
    bp_2_corr = c(-2, 0, 2),
    age_2 = c(-1.5, 0, 1.5),
    age_round = c(-1.50, 0.00, 1.50),
    bmi_2 = c(-2, 0, 2),
    sex = c(0,0,0),
    name_id = c("Id = 1", "Id = 2", "Id = 3"),
    name = c("Participant 1 ", "Participant 2", "Participant 3")
  )
  
  data_fake[, age_effect := predict.gam(age_func, data.frame(age_round = age_round)) ]
  
  # Collect the posterior mean for the different parameters
  b_0 <-inla_spm$summary.fixed[6,1]
  b_bp <-inla_spm$summary.fixed[8,1]
  b_sex <-inla_spm$summary.fixed[7,1]
  b_bmi <-inla_spm$summary.fixed[9,1]
  c <- bri.hyperpar.summary(inla_spm)[3]
  
  #### Residuals ---- 
  #Construct residuals between -1, 1
  
  residuals <- seq(-1.5,1.5, by = 0.01)
  n_residuals <- length(residuals)
  
  resid_data <- data.table(res = rep(residuals, 3),
                           id = c(rep(1, n_residuals), rep(2, n_residuals), rep(3, n_residuals)))
  
  # Merge the risidual data with the fake individuals so every individual gets all residuals
  data_fake_resid<- merge(data_fake, resid_data, on = "id")
  data_fake_resid[, eta_missing := b_0 + b_bp*bp_2_corr + age_effect +b_sex*sex + b_bmi*bmi_2 + c*res]
  data_fake_resid[, p := exp(eta_missing)/(exp(eta_missing)+1)]
  
  # Display the resulting probability of dropping out. 
  q <- ggplot(data = data_fake_resid, aes(x = res, y = p, group = id, colour = as.factor(id), linetype = as.factor(id)))
  q <- q + geom_line()
  q <- q + ylab(TeX("Probability of dropping out")) + xlab(TeX("$\\epsilon_i$")) + theme( legend.position = "bottom")
  q <- q + labs(colour  = "Id", linetype = "Id")
  q <- q + scale_linetype_manual(values = c("solid","dashed", "dotted"))
  q
  
  ggsave("images/assosiation_effect.pdf",
         width = 17,
         height = 10,
         units = "cm")
  
}

