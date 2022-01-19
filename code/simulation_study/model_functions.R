# INLA model functions without config


run_spm <- function(data_sim, c_mu = 0, c_sigma = 1, verbose = F, compute_all = FALSE){
  
  sigma2 = 0.001
  #Prepare the response variables:
  n <- nrow(data_sim)
  y_gaussian <- c(data_sim$bp_3_corr,rep(NA,n))
  m_binomial <- c(rep(NA,n),data_sim$missing)
  joint_response <- list(y_gaussian,m_binomial)
  
  #Must specify number of trials for binomial repsonse variable
  ntrials <- c(rep(NA,n),rep(1,n))
  
  linear_covariates <- data.frame(alpha_0 = c(rep(1,n),rep(NA,n)),
                                  beta_0 = c(rep(NA,n),rep(1,n)),
                                  y_SEX = c(data_sim$sex,rep(NA,n)),
                                  y_BP2 = c(data_sim$bp_2_corr,rep(NA,n)),
                                  y_AGE = c(data_sim$age_2,rep(NA,n)),
                                  y_BMI = c(data_sim$bmi_2,rep(NA,n)),
                                  m_SEX = c(rep(NA,n),data_sim$sex),
                                  m_BP2 = c(rep(NA,n),data_sim$bp_2_corr),
                                  m_AGE = c(rep(NA,n),data_sim$age_2),
                                  m_BMI = c(rep(NA,n),data_sim$bmi_2)
  )
  random_covariates <- data.frame(y_eps1 = c(1:n,rep(NA,n)),
                                  m_eps1 = c(rep(NA,n),1:n))
  
  joint_data <- c(linear_covariates,random_covariates)
  
  joint_data$Y <- joint_response
  
  formula_spm = Y ~ -1 + alpha_0 + y_BP2 + y_AGE + y_SEX + y_BMI +
    beta_0 + m_SEX + m_BP2 +  m_BMI +
    f(m_AGE,model="rw2",constr=T) +
    f(y_eps1,model="iid") +
    f(m_eps1,copy="y_eps1",fixed=F,param=c(c_mu,c_sigma)) #Copy epsilon1 into binomial dropout process with association parameter
  # Use the option param to set mean and precision on Gaussian prior for associaction parameter
  
  link= rep(NA, 2*(n ))
  link[which(is.na(y_gaussian[1:(n )]))] = 1
  link[(n) + which(is.na(m_binomial[((n)+1):(2*(n))]))] = 2
  length(link)
  start_time <- Sys.time()
  if(compute_all == FALSE ){
    joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                       Ntrials=ntrials, verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                       control.predictor = list(link = link))
  }else{
    joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                       Ntrials=ntrials, verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                       control.predictor = list(link = link),
                       control.compute=list(config=TRUE))
    
  }

  
  # Use control.family to fix sigma2. Specified through log-precision.
  
  end_time <- Sys.time()
  return(joint_inla)
}

run_naive_bp <- function(data_sim, verbose = F,compute_all = FALSE){
  n <- nrow(data_sim)
  sigma2 = 0.001 #Fixed to a neglectable value
  
  #Prepare the response variables:
  y_gaussian <- data_sim$bp_3_corr
  
  
  linear_covariates <- list(alpha_0 = rep(1,n),
                            y_SEX = data_sim$sex,
                            y_BP1 = data_sim$bp_2_corr,
                            y_AGE = data_sim$age_2,
                            y_BMI = data_sim$bmi_2,
                            y_eps1 = 1:n)
  naive_data_bp <- linear_covariates
  
  naive_data_bp$Y <- y_gaussian
  
  formula_naive_bp = Y ~ -1 + alpha_0 + y_BP1 + y_AGE + y_SEX + y_BMI +
    f(y_eps1,model="iid") 
  
  start_time <- Sys.time()
  if(compute_all == FALSE){
    joint_inla <- inla(formula_naive_bp, 
                       family="gaussian",
                       data=naive_data_bp, 
                       verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T)))
  }else{
    joint_inla <- inla(formula_naive_bp, 
                       family="gaussian",
                       data=naive_data_bp, 
                       verbose=verbose,
                       control.family=list(list(initial=log(1/sigma2^2),fixed=T)),
                       control.compute=list(config=TRUE))
  }

  return(joint_inla)
}

run_naive_m <- function(data_sim, verbose = F, compute_all = FALSE){
  data<- data_sim
  n <- nrow(data)
  
  linear_covariates_n_m <- list(beta_0 = rep(1,n),
                                m_SEX = data$sex,
                                m_BP1 = data$bp_2_corr,
                                m_AGE = data$age_2,
                                m_BMI = data$bmi_2,
                                m_eps1 =1:n)
  
  naive_data_m <- linear_covariates_n_m
  y_naive_m <- data$missing
  
  naive_data_m$Y <- y_naive_m
  formula_naive_m = Y ~ -1 + beta_0  + m_SEX + m_BP1 + m_BMI + f(m_AGE,model="rw2",constr=T)
  ntrials <- 1
  
  if(compute_all == FALSE){
    inla_naive_m <- inla(formula = formula_naive_m, 
                         family = "binomial", 
                         Ntrials = ntrials, 
                         data = naive_data_m, 
                         verbose=verbose)
  }else{
    inla_naive_m <- inla(formula = formula_naive_m, 
                         family = "binomial", 
                         data = naive_data_m, 
                         Ntrials=ntrials, 
                         verbose=verbose,
                         control.compute=list(config=TRUE))
  }


  
  
  # Use control.family to fix sigma2. Specified through log-precision.
  
  return(inla_naive_m)
}


create_sim_data <- function(spm_param, spm_coverage, naive_param, naive_coverage,  true_values){
  
  spm_long = melt(spm_param, id.vars = c("id"),
                  measure.vars = colnames(spm_param)[which(!colnames(spm_param) %in% c("id"))])
  setnames(spm_long, "value", "value_spm")
  spm_long[true_values, on = "variable", true_value:= value_true]
  
  naive_long = melt(naive_param, id.vars = c("id"),
                    measure.vars = colnames(naive_param)[which(!colnames(naive_param) %in% c("id"))])
  setnames(naive_long, "value", "value_naive")
  
  spm_coverage_long = melt(spm_coverage, id.vars = c("id"),
                           measure.vars = colnames(spm_coverage)[which(!colnames(spm_coverage) %in% c("id"))])
  setnames(spm_coverage_long, "value", "coverage_spm")
  
  naive_coverage_long = melt(naive_coverage, id.vars = c("id"),
                             measure.vars = colnames(naive_coverage)[which(!colnames(naive_coverage) %in% c("id"))])
  setnames(naive_coverage_long, "value", "coverage_naive")
  
  
  simulation_results <- merge(spm_long, naive_long, by = c("id", "variable"), all = TRUE)
  simulation_results <- merge(simulation_results, spm_coverage_long, by = c("id", "variable"), all = TRUE)
  simulation_results <- merge(simulation_results, naive_coverage_long, by = c("id", "variable"), all = TRUE)
  
  simulation_results[, bias_spm := value_spm-true_value]
  simulation_results[, bias_naive := value_naive-true_value]
  
  simulation_results_summary <- simulation_results[, 
                                                   .(true_value = mean(true_value),
                                                     mean_spm = mean(na.omit(value_spm)),
                                                     bias_spm = mean(na.omit(bias_spm)),
                                                     coverage_spm = mean(na.omit(coverage_spm)),
                                                     mean_naive = mean(value_naive),
                                                     bias_naive = mean(bias_naive),
                                                     coverage_naive = mean(coverage_naive) 
                                                   ),
                                                   keyby = .(variable)]
  simulation_results_summary[true_values, on = "variable", true_value := value_true ]
  simulation_results_summary[, diff_bias := abs(bias_naive)-abs(bias_spm)]
  
  retval <- vector(length = 2)
  retval$simulation_results <- simulation_results
  retval$simulation_results_sumary <- simulation_results_summary
  return(retval)
}

plot_sim_results <- function(data, data_summary){
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
  
  q <- ggplot(data, aes(x = value_spm))
  q <- q + geom_density(aes(color = "SPM", linetype = "SPM"))
  q <- q + geom_density( aes(x = value_naive, color = "Naive", linetype = "Naive"))
  q <- q + geom_vline(data = data_summary, aes(xintercept = true_value, color = "True", linetype = "True"))
  q <- q + geom_vline(data = data_summary, aes(xintercept = mean_spm, color = "SPM", linetype = "SPM"))
  q <- q + geom_vline(data = data_summary, aes(xintercept = mean_naive, color = "Naive", linetype = "Naive"))
  q <- q + labs(color  = "Guide name", linetype = "Guide name") + ylab("Density") + xlab(element_blank())
  q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
  q <- q + scale_color_manual(values = c("#F8766D","#00BFC4", "black"))
  return(q + facet_wrap(~variable, ncol = 3, scales = "free", 
                        labeller = labeller(variable = as_labeller(variable_names, label_parsed))))
  
}

q025 <- function(x){
  return(quantile(x, probs = 0.025))
}
q500 <- function(x){
  return(quantile(x, probs = 0.5))
}
q975 <- function(x){
  return(quantile(x, probs = 0.975))
}