# validation func

validation_func_spm <- function(spm, data, n_samples){
  # returns n_sample predictions for SPM and naive model given a dataset 
  
  sigma2 = 0.001
  n_val <- nrow(data)
  
  #### Age effect SPM----
  age_table_spm <- as.data.table(spm$summary.random$m_AGE[, c("ID", "mean")])
  age_table_spm[, age_round := round(ID, 2)]
  age_func_spm <- gam(mean ~s(age_round),data=age_table_spm)
  data[, age_round := round(age_3, 2)]
  data[, age_effect_spm := predict.gam(age_func_spm, data.frame(age_round = age_round)) ]
  
  ############ Samples for parameters from posteriro distribution SPM ###################
  sim_spm <- inla.posterior.sample(n_samples, 
                                   spm, 
                                   selection = list(alpha_0 = 1, 
                                                    y_BP2=1, 
                                                    y_AGE = 1, 
                                                    y_SEX = 1, 
                                                    y_BMI=1, 
                                                    beta_0 = 1, 
                                                    m_SEX=1, 
                                                    m_BP2=1, 
                                                    m_BMI = 1  ))
  
  # Evaluation function
  fun_spm = function(data = NA) {
    eps = rnorm(length(data$bp_3_corr),sd = sqrt(1/theta[2]))
    y_pred <- alpha_0 + y_BP2 * data$bp_3_corr  + y_AGE *data$age_3 + y_SEX*data$sex + y_BMI*data$bmi_3 + eps
    
    eta_m <- beta_0 + m_BP2* data$bp_3_corr + m_SEX*data$sex  +m_BMI*data$bmi_3 + data$age_effect_spm + theta[3]*eps
    m_pred <- exp(eta_m)/(exp(eta_m)+1)
    return (c(y_pred, m_pred))
  }
  print("spm samples ok")
  predictions_spm <- inla.posterior.sample.eval(fun_spm, sim_spm, data = data)
  print("spm ok")
  predictions_spm <- as.data.table(predictions_spm)
  
  predictions_bp_spm <- predictions_spm[1:n_val,]
  predictions_bp_spm[, id := 1:n_val]
  
  
  predictions_m_spm <- predictions_spm[(n_val + 1):(2*n_val),]
  predictions_m_spm[, id := 1:n_val]
  
  
  validation_data <- vector(mode= "list")
  validation_data$predictions_bp <-predictions_bp_spm
  validation_data$predictions_m <- predictions_m_spm
  return(validation_data)
}

validation_func_naive <- function(naive_bp, naive_m, data, n_samples){
  # returns n_sample predictions for SPM and naive model given a dataset 
  
  sigma2 = 0.001
  n_val <- nrow(data)
  
  #### Age effect naive ----
  age_table_n <- as.data.table(naive_m$summary.random$m_AGE[, c("ID", "mean")])
  age_table_n[, age_round := round(ID, 2)]
  
  age_func_n <- gam(mean ~s(age_round),data=age_table_n)
  data[, age_round := round(age_3, 2)]
  data[, age_effect_n := predict.gam(age_func_n, data.frame(age_round = age_round)) ]
  
  ##### Samples for parameters from posteriro distribution naive ----
  
  sim_n_bp <- inla.posterior.sample(n_samples, naive_bp, selection = list(alpha_0 = 1, y_BP1=1, y_AGE = 1, y_SEX = 1, y_BMI=1  ))
  sim_n_m <- inla.posterior.sample(n_samples, naive_m, selection = list( beta_0 = 1, m_SEX=1, m_BP1=1, m_BMI = 1  ))
  print("naive sample ok")
  # Evaluation function
  fun_n_bp = function(data = NA) {
    eps = rnorm(length(data$bp_3_corr),sd = sqrt(1/theta[1]))
    y_pred <- alpha_0 + y_BP1 * data$bp_3_corr  + y_AGE *data$age_3 + y_SEX*data$sex + y_BMI*data$bmi_3 + eps
    
    return (y_pred)
  }
  
  fun_n_m = function(data = NA) {
    eta_m <- beta_0 + m_BP1* data$bp_3_corr + m_SEX*data$sex  +m_BMI*data$bmi_3 + data$age_effect_n
    m_pred <- exp(eta_m)/(exp(eta_m)+1)
    return (c(m_pred))
  }
  
  
  predictions_bp_n = inla.posterior.sample.eval(fun_n_bp, sim_n_bp, data = data)
  predictions_bp_n<- as.data.table(predictions_bp_n)
  predictions_bp_n[, id := 1:n_val]
  
  predictions_m_n = inla.posterior.sample.eval(fun_n_m, sim_n_m, data = data)
  predictions_m_n<- as.data.table(predictions_m_n)
  predictions_m_n[, id := 1:n_val]
  print("naive ok" )
  validation_data <- vector(mode= "list")
  validation_data$predictions_bp <- predictions_bp_n
  validation_data$predictions_m <- predictions_m_n
  return(validation_data)
}

validation_func <- function(spm, naive_bp, naive_m, data, n_samples){
  # returns n_sample predictions for SPM and naive model given a dataset 
  
  sigma2 = 0.001
  n_val <- nrow(data)

  #### Age effect SPM----
  age_table_spm <- as.data.table(spm$summary.random$m_AGE[, c("ID", "mean")])
  age_table_spm[, age_round := round(ID, 2)]
  age_func_spm <- gam(mean ~s(age_round),data=age_table_spm)
  data[, age_round := round(age_3, 2)]
  data[, age_effect_spm := predict.gam(age_func_spm, data.frame(age_round = age_round)) ]
  
  ############ Samples for parameters from posteriro distribution SPM ###################
  sim_spm <- inla.posterior.sample(n_samples, 
                                   spm, 
                                   selection = list(alpha_0 = 1, 
                                                    y_BP2=1, 
                                                    y_AGE = 1, 
                                                    y_SEX = 1, 
                                                    y_BMI=1, 
                                                    beta_0 = 1, 
                                                    m_SEX=1, 
                                                    m_BP2=1, 
                                                    m_BMI = 1  ))
  
  # Evaluation function
  fun_spm = function(data = NA) {
    eps = rnorm(length(data$bp_3_corr),sd = sqrt(1/theta[2]))
    y_pred <- alpha_0 + y_BP2 * data$bp_3_corr  + y_AGE *data$age_3 + y_SEX*data$sex + y_BMI*data$bmi_3 + eps

    eta_m <- beta_0 + m_BP2* data$bp_3_corr + m_SEX*data$sex  +m_BMI*data$bmi_3 + data$age_effect_spm + theta[3]*eps
    m_pred <- exp(eta_m)/(exp(eta_m)+1)
    return (c(y_pred, m_pred))
  }
  print("spm samples ok")
  predictions_spm <- inla.posterior.sample.eval(fun_spm, sim_spm, data = data)
  print("spm ok")
  predictions_spm <- as.data.table(predictions_spm)

  predictions_bp_spm <- predictions_spm[1:n_val,]
  predictions_bp_spm[, id := 1:n_val]
 

  predictions_m_spm <- predictions_spm[(n_val + 1):(2*n_val),]
  predictions_m_spm[, id := 1:n_val]
 

  #### Age effect naive ----
  age_table_n <- as.data.table(naive_m$summary.random$m_AGE[, c("ID", "mean")])
  age_table_n[, age_round := round(ID, 2)]
  
  age_func_n <- gam(mean ~s(age_round),data=age_table_n)
  data[, age_round := round(age_3, 2)]
  data[, age_effect_n := predict.gam(age_func_n, data.frame(age_round = age_round)) ]
  
  ##### Samples for parameters from posteriro distribution naive ----
  
  sim_n_bp <- inla.posterior.sample(n_samples, naive_bp, selection = list(alpha_0 = 1, y_BP1=1, y_AGE = 1, y_SEX = 1, y_BMI=1  ))
  sim_n_m <- inla.posterior.sample(n_samples, naive_m, selection = list( beta_0 = 1, m_SEX=1, m_BP1=1, m_BMI = 1  ))
  print("naive sample ok")
  # Evaluation function
  fun_n_bp = function(data = NA) {
    eps = rnorm(length(data$bp_3_corr),sd = sqrt(1/theta[1]))
    y_pred <- alpha_0 + y_BP1 * data$bp_3_corr  + y_AGE *data$age_3 + y_SEX*data$sex + y_BMI*data$bmi_3 + eps

    return (y_pred)
  }

  fun_n_m = function(data = NA) {
    eta_m <- beta_0 + m_BP1* data$bp_3_corr + m_SEX*data$sex  +m_BMI*data$bmi_3 + data$age_effect_n
    m_pred <- exp(eta_m)/(exp(eta_m)+1)
    return (c(m_pred))
  }

  
  predictions_bp_n = inla.posterior.sample.eval(fun_n_bp, sim_n_bp, data = data)
  predictions_bp_n<- as.data.table(predictions_bp_n)
  predictions_bp_n[, id := 1:n_val]

  predictions_m_n = inla.posterior.sample.eval(fun_n_m, sim_n_m, data = data)
  predictions_m_n<- as.data.table(predictions_m_n)
  predictions_m_n[, id := 1:n_val]
  print("naive ok" )
  validation_data <- vector(mode= "list")
  validation_data$predictions_bp_spm <-predictions_bp_spm
  validation_data$predictions_m_spm <- predictions_m_spm
  validation_data$predictions_bp_n <- predictions_bp_n
  validation_data$predictions_m_n <- predictions_m_n
  return(validation_data)
}

crps_func_general <- function(predictions_bp, data_val, n_samples){
  # Returns CRPS score for spm and naive model
  predictions_bp[data_val, on = ("id"), bp_corr:= bp_4_corr]
  predictions_bp[data_val, on = ("id"), missing:= missing]
  predictions_bp_compleat <- na.omit(predictions_bp)
  
  crps <- crps_sample(as.vector(predictions_bp_compleat$bp_corr), 
                          dat = as.matrix(predictions_bp_compleat[, 1:(n_samples)]))
  
  crps_data <- data.table(id = predictions_bp_compleat$id, 
                          bp_corr = predictions_bp_compleat$bp_corr,
                          missing = predictions_bp_compleat$missing, 
                          crps = crps)
  crps_total <- copy(crps_data)
  crps_total[, missing := -1]
  
  crps <- as.data.table(rbindlist(list(crps_data, crps_total)))
  
  crps_agg <- crps[, .(crps = mean(crps)
  ),
  keyby = .(missing)]
  return(crps_agg)
}


crps_func <- function(predictions_bp_spm, predictions_bp_n, data_val, n_samples){
  # Returns CRPS score for spm and naive model
  predictions_bp_spm[data_val, on = ("id"), bp_corr:= bp_4_corr]
  predictions_bp_spm[data_val, on = ("id"), missing:= missing]
  predictions_bp_spm_compleat <- na.omit(predictions_bp_spm)
  
  crps_spm <- crps_sample(as.vector(predictions_bp_spm_compleat$bp_corr), 
                          dat = as.matrix(predictions_bp_spm_compleat[, 1:(n_samples)]))
  
  predictions_bp_n[data_val, on = ("id"), bp_corr:= bp_4_corr]
  predictions_bp_n[data_val, on = ("id"), missing:= missing]
  predictions_bp_n_compleat <- na.omit(predictions_bp_n)
  
  crps_n <- crps_sample(as.vector(predictions_bp_n_compleat$bp_corr), 
                        dat = as.matrix(predictions_bp_n_compleat[,1:(n_samples)]))
  crps_data <- data.table(id = predictions_bp_spm_compleat$id, 
                          bp_corr = predictions_bp_spm_compleat$bp_corr,
                          missing = predictions_bp_spm_compleat$missing, 
                          crps_spm = crps_spm, 
                          crps_n = crps_n)
  crps_total <- copy(crps_data)
  crps_total[, missing := -1]
  
  crps <- as.data.table(rbindlist(list(crps_data, crps_total)))
  
  crps_agg <- crps[, .(crps_spm = mean(crps_spm), 
                       crps_n = mean(crps_n)
                       ),
                   keyby = .(missing)]
  return(crps_agg)
  
  
}

brier_func_general <- function(predictions_m, data_val){
  #returns brier score for spm and naive model
  predictions_m_long <- melt(predictions_m, id.vars= c("id"))
  setnames(predictions_m_long, c("value", "variable"), c("m_pred", "row_id"))
  
  predictions_m_agg <- predictions_m_long[, .(m_pred = mean(m_pred)),
                                                  keyby = .(id)]
  
  predictions_m_agg[data_val, on = "id", missing := missing]
  
  
  predictions_m <- predictions_m_agg
  
  predictions_m[, diff_squared := (m_pred-missing)^2]
  
  predictions_m_total <- copy(predictions_m)
  predictions_m_total[, missing := -1]
  
  predictions_m <- rbindlist(list(predictions_m, predictions_m_total))
  predictions_m_agg<- predictions_m[, .(brier = mean(diff_squared)), 
                                    keyby = .(missing)]
  return(predictions_m_agg)
}

brier_func <- function(predictions_m_spm, predictions_m_n, data_val){
  #returns brier score for spm and naive model
  predictions_m_long_spm <- melt(predictions_m_spm, id.vars= c("id"))
  setnames(predictions_m_long_spm, c("value", "variable"), c("m_pred", "row_id"))
  predictions_m_long_n <- melt(predictions_m_n, id.vars= c("id"))
  setnames(predictions_m_long_n, c("value", "variable"), c("m_pred", "row_id"))
  
  predictions_m_spm_agg <- predictions_m_long_spm[, .(m_pred_spm = mean(m_pred)),
                                                  keyby = .(id)]
  predictions_m_n_agg <- predictions_m_long_n[, .(m_pred_n = mean(m_pred)),
                                              keyby = .(id)]
  
  predictions_m_spm_agg[data_val, on = "id", missing := missing]
  

  predictions_m <- predictions_m_spm_agg
  predictions_m[predictions_m_n_agg, on = "id", m_pred_n := m_pred_n]
  
  predictions_m[, diff_squared_spm := (m_pred_spm-missing)^2]
  predictions_m[, diff_squared_n := (m_pred_n-missing)^2]
  
  predictions_m_total <- copy(predictions_m)
  predictions_m_total[, missing := -1]
  
  predictions_m <- rbindlist(list(predictions_m, predictions_m_total))
  predictions_m_agg<- predictions_m[, .(brier_spm = mean(diff_squared_spm),
                    brier_n = mean(diff_squared_n)), 
                        keyby = .(missing)]
  return(predictions_m_agg)
}



pred_spm_given_m <- function(data_sim, true_values, n){
  sigma2<- 0.001
  a_0  = true_values[1,2][[1]]
  a_bp = true_values[2,2][[1]]
  a_age = true_values[3,2][[1]]
  a_sex = true_values[4,2][[1]]
  a_bmi = true_values[5,2][[1]]
  b_0 = true_values[6,2][[1]]
  b_sex = true_values[7,2][[1]]
  b_bp = true_values[8,2][[1]]
  b_bmi = true_values[9,2][[1]]
  
  sigma_age <- true_values[10,2][[1]]
  sigma <- true_values[11,2][[1]]
  c <- true_values[12,2][[1]]
  
  
  precision_age_mean <-true_values[13,2][[1]]
  precision_age_sd <- true_values[14,2][[1]]
  precision_eps_mean <- true_values[15,2][[1]]
  precision_eps_sd <- true_values[16,2][[1]]
  
  mean_c <- true_values[17,2][[1]]
  sd_c <- true_values[18,2][[1]]
  
  
  #### SPM predict BP_f given m #### 
  sigma2 = 0.001
  y_gaussian <- c(rep(NA,n),rep(NA,n)) #set y to na
  m_binomial <- c(rep(NA,n),data_sim$missing)
  joint_response <- list(y_gaussian,m_binomial)
  
  ntrials <- c(rep(NA,n),rep(1,n))
  
  linear_covariates <- data.frame(alpha_0 = c(rep(1,n),rep(NA,n)),
                                  beta_0 = c(rep(NA,n),data_sim$intercept),
                                  y_SEX = c(data_sim$sex,rep(NA,n)),
                                  y_BP2 = c(data_sim$bp_3_corr,rep(NA,n)),
                                  y_AGE = c(data_sim$age_3,rep(NA,n)),
                                  y_BMI = c(data_sim$bmi_3,rep(NA,n)),
                                  m_SEX = c(rep(NA,n),data_sim$sex),
                                  m_BP2 = c(rep(NA,n),data_sim$bp_3_corr),
                                  m_AGE = c(rep(NA,n),data_sim$age_3),
                                  m_BMI = c(rep(NA,n),data_sim$bmi_3)
  )
  random_covariates <- data.frame(y_eps1 = c(1:n,rep(NA,n)),
                                  m_eps1 = c(rep(NA,n),1:n))
  
  joint_data <- c(linear_covariates,random_covariates)
  
  joint_data$Y <- joint_response
  
  #prec.prior_age <- list(prec = list(prior = "gaussian", param = c(precision_age_mean, precision_age_sd), fixed = TRUE))
  prec.prior_eps <- list(prec = list(prior = "gaussian", param = c(precision_eps_mean, precision_eps_sd), fixed = TRUE))
  formula_spm = Y ~ -1 + alpha_0 + y_BP2 + y_AGE + y_SEX + y_BMI +
    beta_0 + m_SEX + m_BP2 +  m_BMI + 
    f(y_eps1,model="iid", hyper = prec.prior_eps) +
    f(m_eps1,copy="y_eps1",fixed=T,param=c(mean_c,sd_c)) 
  
  precision = 1000000
  
  link= rep(NA, 2*(n ))
  link[which(is.na(y_gaussian[1:(n )]))] = 1
  link[(n) + which(is.na(m_binomial[((n)+1):(2*(n))]))] = 2
  length(link)
  start_time <- Sys.time()
  
  # Specify very narrow prior distributions for parameters
  joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                     Ntrials=ntrials, verbose=F,
                     control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                     control.predictor = list(link = link), 
                     control.fixed=list(mean.intercept = 1, prec.intercept = precision, 
                                        mean=list(alpha_0= a_0 ,y_BP2= a_bp,  y_AGE= a_age, y_SEX= a_sex,  y_BMI= a_bmi, 
                                                  beta_0= 1,  m_SEX= b_sex,  m_BP2= b_bp,   m_BMI=b_bmi , default=0), 
                                        prec = list(alpha_0= precision ,y_BP2= precision,  y_AGE= precision, y_SEX= precision,  y_BMI=precision , 
                                                    beta_0= precision,  m_SEX=precision ,  m_BP2=precision,   m_BMI=precision , default=0.01))
  )
  return(joint_inla)
}

pred_spm_bp <- function(sim_data, true_values,n){
  sigma2<- 0.001
  a_0  = true_values[1,2][[1]]
  a_bp = true_values[2,2][[1]]
  a_age = true_values[3,2][[1]]
  a_sex = true_values[4,2][[1]]
  a_bmi = true_values[5,2][[1]]
  b_0 = true_values[6,2][[1]]
  b_sex = true_values[7,2][[1]]
  b_bp = true_values[8,2][[1]]
  b_bmi = true_values[9,2][[1]]
  
  sigma_age <- true_values[10,2][[1]]
  sigma <- true_values[11,2][[1]]
  c <- true_values[12,2][[1]]
  
  
  precision_age_mean <-true_values[13,2][[1]]
  precision_age_sd <- true_values[14,2][[1]]
  precision_eps_mean <- true_values[15,2][[1]]
  precision_eps_sd <- true_values[16,2][[1]]
  
  mean_c <- true_values[17,2][[1]]
  sd_c <- true_values[18,2][[1]]
  # spm predict BP without m
  
  ntrials <- c(rep(NA,n),rep(1,n))
  
  linear_covariates <- data.frame(alpha_0 = c(rep(1,n),rep(NA,n)),
                                  beta_0 = c(rep(NA,n),data_sim$intercept),
                                  y_SEX = c(data_sim$sex,rep(NA,n)),
                                  y_BP2 = c(data_sim$bp_3_corr,rep(NA,n)),
                                  y_AGE = c(data_sim$age_3,rep(NA,n)),
                                  y_BMI = c(data_sim$bmi_3,rep(NA,n)),
                                  m_SEX = c(rep(NA,n),data_sim$sex),
                                  m_BP2 = c(rep(NA,n),data_sim$bp_3_corr),
                                  m_AGE = c(rep(NA,n),data_sim$age_3),
                                  m_BMI = c(rep(NA,n),data_sim$bmi_3)
  )
  random_covariates <- data.frame(y_eps1 = c(1:n,rep(NA,n)),
                                  m_eps1 = c(rep(NA,n),1:n))
  
  joint_data_n <- c(linear_covariates,random_covariates)
  y_gaussian_n <- c(rep(NA,n),rep(NA,n)) #set y to na
  m_binomial_n <- c(rep(NA,n),rep(NA,n))
  joint_response_n <- list(y_gaussian_n,m_binomial_n)
  joint_data_n$Y <- joint_response_n
  
  prec.prior_eps <- list(prec = list(prior = "gaussian", param = c(precision_eps_mean, precision_eps_sd), fixed = TRUE))
  formula_spm = Y ~ -1 + alpha_0 + y_BP2 + y_AGE + y_SEX + y_BMI +
    beta_0 + m_SEX + m_BP2 +  m_BMI + 
    f(y_eps1,model="iid", hyper = prec.prior_eps) +
    f(m_eps1,copy="y_eps1",fixed=T,param=c(mean_c,sd_c)) 
  
  precision = 1000000
  
  link= rep(NA, 2*(n ))
  link[which(is.na(y_gaussian_n[1:(n )]))] = 1
  link[(n) + which(is.na(m_binomial_n[((n)+1):(2*(n))]))] = 2
  length(link)
  
  joint_inla_n <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data_n,
                       Ntrials=ntrials, verbose=F,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                       control.predictor = list(link = link), 
                       control.fixed=list(mean.intercept = 1, prec.intercept = precision, 
                                          mean=list(alpha_0= a_0 ,y_BP2= a_bp,  y_AGE= a_age, y_SEX= a_sex,  y_BMI= a_bmi, 
                                                    beta_0= 1,  m_SEX= b_sex,  m_BP2= b_bp,   m_BMI=b_bmi , default=0), 
                                          prec = list(alpha_0= precision ,y_BP2= precision,  y_AGE= precision, y_SEX= precision,  y_BMI=precision , 
                                                      beta_0= precision,  m_SEX=precision ,  m_BP2=precision,   m_BMI=precision , default=0.01))
  )
  
  summary(joint_inla_n)
  true_values
  return(joint_inla_n)
}






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
# #### Support functions ----
# q025 <- function(x){
#   return(quantile(x, probs = 0.025))
# }
# q500 <- function(x){
#   return(quantile(x, probs = 0.5))
# }
# q975 <- function(x){
#   return(quantile(x, probs = 0.975))
# }
# 
# ############  SPM ###################
# #### Read the model and data ----
# 
# spm <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT23.RData")
# naive_bp <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/inla_naive_bp_true_missing.RData" )
# naive_m <- readRDS(file =  "/home/aurorach/data_HUNT_aurora/master/inla_naive_m_true_missing.RData")
# 
# data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")
# data <- data[1:10,]
# n_samples <- 3 #1000 #300