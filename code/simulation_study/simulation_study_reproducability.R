library(INLA)
library(brinla)
#INLA:::inla.dynload.workaround()
library(tidyverse)
library(data.table)
library(latex2exp)
library(tikzDevice)
library(lubridate)
library(glue)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(gridtext)


source("code/simulation_study/model_functions.R") 



simulation_study_reproducability_func <- function(data, n_sim = 100, which_parameters = "SPM", is_mnar = TRUE, dir){
  n = nrow(data)
  
  if(which_parameters == "SPM"){
    # Read parameters inla spm
    inla_spm<-readRDS(file = glue::glue("{dir}/spm_inla_HUNT2.RData"))
    a_0  <- inla_spm$summary.fixed[[1]][1]
    a_bp <- inla_spm$summary.fixed[[1]][2]
    a_age <- inla_spm$summary.fixed[[1]][3]
    a_sex <- inla_spm$summary.fixed[[1]][4]
    a_bmi <- inla_spm$summary.fixed[[1]][5]
    b_0 <- inla_spm$summary.fixed[[1]][6]
    b_sex <- inla_spm$summary.fixed[[1]][7]
    b_bp <- inla_spm$summary.fixed[[1]][8]
    b_bmi <- inla_spm$summary.fixed[[1]][9]
    
    sigma_age <- bri.hyperpar.summary(inla_spm)[1]
    sigma <- bri.hyperpar.summary(inla_spm)[2]
    if(is_mnar){
      c <-  bri.hyperpar.summary(inla_spm)[3]
    }else{
      c <- 0
    }
    summary_fixed <- inla_spm$summary.fixed[[1]]
    summary_hyper <- bri.hyperpar.summary(inla_spm)[1:3]
    summary_hyper[3]<- c
    
    age_table <- as.data.table(inla_spm$summary.random$m_AGE[, c("ID", "mean", "0.025quant", "0.975quant")])
    age_table[, age_round := round(ID, 2)]
    variable = c("a_0", "a_bp", "a_age",  "a_sex", "a_bmi", "b_0", "b_sex", "b_bp", "b_bmi", "sigma_age", "sigma_eps", "c")
    true_values <- data.table(variable = variable,
                              value_true = c(inla_spm$summary.fixed[[1]], bri.hyperpar.summary(inla_spm)[1:2], c))
    
    rm(inla_spm)
  }else if(which_parameters == "naive"){
    #Read parameters inla naive models
    inla_naive_bp<-readRDS(file = glue::glue("{dir}/naive_bp_inla_HUNT2.RData"))
    inla_naive_m<-readRDS(file = glue::glue("{dir}/naive_m_inla_HUNT2.RData"))
    
    a_0  <- inla_naive_bp$summary.fixed[[1]][1]
    a_bp <- inla_naive_bp$summary.fixed[[1]][2]
    a_age <- inla_naive_bp$summary.fixed[[1]][3]
    a_sex <- inla_naive_bp$summary.fixed[[1]][4]
    a_bmi <- inla_naive_bp$summary.fixed[[1]][5]
    b_0 <- inla_naive_m$summary.fixed[[1]][1]
    b_sex <- inla_naive_m$summary.fixed[[1]][2]
    b_bp <- inla_naive_m$summary.fixed[[1]][3]
    b_bmi <- inla_naive_m$summary.fixed[[1]][4]
    
    
    sigma_age <- bri.hyperpar.summary(inla_naive_m)[1]
    sigma <- bri.hyperpar.summary(inla_naive_bp)[1]
    if(is_mnar){
      c <- 0.705
    }else{
      c <- 0
    }
    
    # Create list of variable names
    variable = c("a_0", "a_bp", "a_age",  "a_sex", "a_bmi", "b_0", "b_sex", "b_bp", "b_bmi", "sigma_age", "sigma_eps", "c")
    
    # Store the true values of the parameters
    true_values <- data.table(variable = variable,
                              value_true = c(inla_naive_bp$summary.fixed[[1]],
                                             inla_naive_m$summary.fixed[[1]], 
                                             bri.hyperpar.summary(inla_naive_m)[1], 
                                             bri.hyperpar.summary(inla_naive_bp)[1],
                                             c))
    
    # Store summaries of the true model parameters
    summary_fixed <- c(inla_naive_bp$summary.fixed[[1]], 
                       inla_naive_m$summary.fixed[[1]])
    
    summary_hyper <- c(bri.hyperpar.summary(inla_naive_m)[1], 
                       bri.hyperpar.summary(inla_naive_bp)[1],
                       c)
    #Get age effect
    age_table <- as.data.table(inla_naive_m$summary.random$m_AGE[, c("ID", "mean", "0.025quant", "0.975quant")])
    age_table[, age_round := round(ID, 2)]
    
    
  }
  
  # Directory for storing results ----
  today <- today()
  if(!dir.exists(glue("results/simulations/{today}"))){
    dir.create(glue("results/simulations/{today}"))
  }
  
  if(is_mnar){
    name <- paste0(which_parameters, "_parameters_MNAR")
  }else{
    name <- paste0(which_parameters, "_parameters_MAR")
  }
  
  if(!dir.exists(glue("results/simulations/{today}/{name}"))){
    dir.create(glue("results/simulations/{today}/{name}"))
  }
  write.csv(as.data.table(true_values), file = glue("results/simulations/{today}/{name}/true_values"), row.names = FALSE )
  
  # Tables for simulation results ---- 
  param_est_spm <- data.table(
    a_0 = rep(0, n_sim),
    a_bp = rep(0, n_sim),
    a_age = rep(0, n_sim),
    a_sex = rep(0, n_sim),
    a_bmi = rep(0, n_sim),
    b_0 = rep(0, n_sim),
    b_sex = rep(0, n_sim),
    b_bp = rep(0, n_sim),
    b_bmi = rep(0, n_sim),
    sigma_age = rep(0, n_sim),
    sigma_eps = rep(0, n_sim),
    c = rep(0, n_sim),
    id = 1:n_sim
  )
  
  coverage_spm <- data.table(
    a_0 = rep(-1, n_sim),
    a_bp = rep(-1, n_sim),
    a_age = rep(-1, n_sim),
    a_sex = rep(-1, n_sim),
    a_bmi = rep(-1, n_sim),
    b_0 = rep(-1, n_sim),
    b_sex = rep(-1, n_sim),
    b_bp = rep(-1, n_sim),
    b_bmi = rep(-1, n_sim),
    sigma_age = rep(-1, n_sim),
    sigma_eps = rep(-1, n_sim),
    c = rep(-1, n_sim),
    id = 1:n_sim
  )
  
  param_est_naive <- data.table(
    a_0 = rep(0, n_sim),
    a_bp = rep(0, n_sim),
    a_age = rep(0, n_sim),
    a_sex = rep(0, n_sim),
    a_bmi = rep(0, n_sim),
    b_0 = rep(0, n_sim),
    b_sex = rep(0, n_sim),
    b_bp = rep(0, n_sim),
    b_bmi = rep(0, n_sim),
    sigma_age = rep(0, n_sim),
    sigma_eps = rep(0, n_sim),
    id = 1:n_sim
  )
  
  coverage_naive <- data.table(
    a_0 = rep(-1, n_sim),
    a_bp = rep(-1, n_sim),
    a_age = rep(-1, n_sim),
    a_sex = rep(-1, n_sim),
    a_bmi = rep(-1, n_sim),
    b_0 = rep(-1, n_sim),
    b_sex = rep(-1, n_sim),
    b_bp = rep(-1, n_sim),
    b_bmi = rep(-1, n_sim),
    sigma_age = rep(-1, n_sim),
    sigma_eps = rep(-1, n_sim),
    id = 1:n_sim
  )
  
  # Time use and failiours  -----
  time_use_spm <- rep(0, n_sim)
  time_use_naive_bp <- rep(0, n_sim)
  time_use_naive_m <- rep(0, n_sim)
  nr_failes_spm <- 0
  nr_failes_hyper <- 0
  sigma2 = 0.001
  
  # Store explanatory variabels ----
  expl_variables <- colnames(data)[which(!colnames(data) %in% c("bp_3_corr", "id", "missing"))]
  
  data_exp <- copy(data[, ..expl_variables])
  data_exp[, age_round:= round(age_2, 2)]
  data_exp[age_table, on="age_round", age_effekt:= mean]
  
  # Start simulations ----
  
  for(i in 1:n_sim){
    
    start_overall<- Sys.time()
    
    #create data sim
    data_sim <- data_exp
    data_sim[, y_eps := rnorm(.N, 0, sigma)]
    data_sim[, y_eps_2 := rnorm(.N, 0, sigma2)]
    data_sim[, eps:= y_eps + y_eps_2]
    data_sim[, bp_3_corr := a_0 + a_bp*bp_2_corr + a_age*age_2 + a_sex*sex + a_bmi*bmi_2 + eps]
    
    data_sim[, eta_p := b_0 + b_bp*bp_2_corr + age_effekt + b_sex*sex + b_bp*bmi_2 + c*y_eps]
    data_sim[, p_i := exp(eta_p)/(exp(eta_p)+1)]
    data_sim[, missing := rbinom(.N, 1, p_i)]
    
    data_sim[missing == 1, bp_3_corr:= NA]
    print(glue::glue("Run {i}"))
    
    # SPM 
    start_time <- Sys.time()
    spm <- try(run_spm(data_sim))
    end_time <- Sys.time()
    if (class(spm) == "try-error"){
      nr_failes_spm = nr_failes_spm + 1
    } else{
      try(print(spm$summary.fixed))
      no_hyper <- try(print(bri.hyperpar.summary(spm)))
      #store parameter estimates
      
      for (j in 1:9){
        param_est_spm[i, j] <- spm$summary.fixed[[1]][j]
      }
      #check coverage
      for (j  in 1:9){
        coverage_spm[i, j] <- as.integer(summary_fixed[j]>=spm$summary.fixed[[3]][j]
                                         & summary_fixed[j] <=spm$summary.fixed[[5]][j])
      }
      
      if((!class(no_hyper) == "try-error")[1]){
        print(class(no_hyper))
        for (j in 1:3){
          param_est_spm[i, j+9] <- bri.hyperpar.summary(spm)[j]
        }
        
        for (j in 1:3){
          coverage_spm[i, j+9] <- as.integer(summary_hyper[j] >= bri.hyperpar.summary(spm)[j +3*2]
                                             & summary_hyper[j] <= bri.hyperpar.summary(spm)[j +3*4])
        }
      }else{
        nr_failes_hyper <- nr_failes_hyper + 1
        for (j in 1:3){
          param_est_spm[i, j+9] <- NA
        }
        
        for (j in 1:3){
          coverage_spm[i, j+9] <-NA
        }
      }
      
    }
    time_use_spm[i] <- end_time - start_time
    print(glue::glue("Time spm round {i}: {end_time - start_time}, nr failes_spm: {nr_failes_spm}, nr_failes_hyper: {nr_failes_hyper}"))
    
    
    
    # NAIVE 
    start_time <- Sys.time()
    naive_bp <- run_naive_bp(data_sim)
    end_time <- Sys.time()
    time_use_naive_bp[i] <- end_time - start_time
    print(glue::glue("Time naive bp round {i}: {end_time - start_time}"))
    #store parameter estimates
    for (j in 1:5){
      param_est_naive[i, j] <- naive_bp$summary.fixed[[1]][j]
    }
    
    param_est_naive[i, 11] <- bri.hyperpar.summary(naive_bp)[1]
    
    #check coverage
    for (j  in 1:5){
      coverage_naive[i, j] <- as.integer(summary_fixed[j]>=naive_bp$summary.fixed[[3]][j]
                                         & summary_fixed[j] <=naive_bp$summary.fixed[[5]][j])
    }
    
    coverage_naive[i, 11] <- as.integer(summary_hyper[2] >= bri.hyperpar.summary(naive_bp)[3]
                                        & summary_hyper[2] <= bri.hyperpar.summary(naive_bp)[5])
    start_time <- Sys.time()
    naive_missing <- run_naive_m(data_sim)
    end_time <- Sys.time()
    time_use_naive_m[i] <- end_time - start_time
    print(glue::glue("Time naive missing round {i}: {end_time - start_time}"))
    #store parameter estimates
    for (j in 6:9){
      param_est_naive[i, j] <- naive_missing$summary.fixed[[1]][j-5]
    }
    
    param_est_naive[i, 10] <- bri.hyperpar.summary(naive_missing)[1]
    
    #check coverage
    for (j  in 6:9){
      coverage_naive[i, j] <- as.integer(summary_fixed[j]>=naive_missing$summary.fixed[[3]][j-5]
                                         & summary_fixed[j] <=naive_missing$summary.fixed[[5]][j-5])
    }
    
    coverage_naive[i, 10] <- as.integer(summary_hyper[1] >= bri.hyperpar.summary(naive_missing)[3]
                                        & summary_hyper[1] <= bri.hyperpar.summary(naive_missing)[5])
    
    end_overall<- Sys.time()
    
    print(glue::glue("Time entire round {i}: {end_overall - start_overall}"))
    if(!i %% 50){
      write.csv(param_est_spm, file = glue("results/simulations/{today}/{name}/{i}_parameter_est_spm"), row.names = FALSE)
      write.csv(param_est_naive, file = glue("results/simulations/{today}/{name}/{i}_parameter_est_naive"), row.names = FALSE)
      write.csv(coverage_spm, file = glue("results/simulations/{today}/{name}/{i}_coverage_spm"), row.names = FALSE)
      write.csv(coverage_naive, file =glue("results/simulations/{today}/{name}/{i}_coverage_naive"), row.names = FALSE)
      write.csv(time_use_spm, file = glue("results/simulations/{today}/{name}/{i}_time_use_spm"), row.names = FALSE)
      write.csv(time_use_naive_bp, file = glue("results/simulations/{today}/{name}/{i}_time_use_bp"), row.names = FALSE)
      write.csv(time_use_naive_m, file = glue("results/simulations/{today}/{name}/{i}_time_use_ms"), row.names = FALSE)
    }
  }
  
  write.csv(param_est_spm, file = glue("results/simulations/{today}/{name}/{n_sim}_parameter_est_spm"), row.names = FALSE)
  write.csv(param_est_naive, file = glue("results/simulations/{today}/{name}/{n_sim}_parameter_est_naive"), row.names = FALSE)
  write.csv(coverage_spm, file = glue("results/simulations/{today}/{name}/{n_sim}_coverage_spm"), row.names = FALSE)
  write.csv(coverage_naive, file =glue( "results/simulations/{today}/{name}/{n_sim}_coverage_naive"), row.names = FALSE)
  write.csv(time_use_spm, file = glue("results/simulations/{today}/{name}/{n_sim}_time_use_spm"), row.names = FALSE)
  write.csv(time_use_naive_bp, file = glue("results/simulations/{today}/{name}/{n_sim}_time_use_bp"), row.names = FALSE)
  write.csv(time_use_naive_m, file = glue("results/simulations/{today}/{name}/{n_sim}_time_use_m"), row.names = FALSE)
  
}

### Functions ----

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

# Read data ----
data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")
#data <- data[1:1000,]
# number of simulations
n_sim = 100

# Run simulation studies ----
simulation_study_reproducability_func(data, n_sim, which_parameters = "SPM", is_mnar = TRUE, dir = "/home/aurorach/data_HUNT_aurora/master")

simulation_study_reproducability_func(data, n_sim, which_parameters = "SPM", is_mnar = FALSE, dir = "/home/aurorach/data_HUNT_aurora/master")

simulation_study_reproducability_func(data, n_sim, which_parameters= "naive", is_mnar = TRUE, dir = "/home/aurorach/data_HUNT_aurora/master")

simulation_study_reproducability_func(data, n_sim, which_parameters = "naive", is_mnar = FALSE, dir = "/home/aurorach/data_HUNT_aurora/master")


# Read simulation study results ----
dir <- "results/simulations/2021-12-17"
### SPM estimates MNAR----

spm_est_coverage_spm_mnar<-  as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MNAR/100_coverage_spm")))
spm_est_coverage_naive_mnar<-  as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MNAR/100_coverage_naive")))

spm_est_parameter_est_spm_mnar <- as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MNAR/100_parameter_est_spm")))
spm_est_parameter_est_naive_mnar <-  as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MNAR/100_parameter_est_naive")))

spm_est_true_values_mnar <- as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MNAR/true_values")))


spm_est_mnar <- create_sim_data(spm_est_parameter_est_spm_mnar, 
                                spm_est_coverage_spm_mnar, 
                                spm_est_parameter_est_naive_mnar, 
                                spm_est_coverage_naive_mnar,
                                spm_est_true_values_mnar)

spm_est_mnar$simulation_results_sumary

q_spm_est_mnar<- plot_sim_results(spm_est_mnar$simulation_results, spm_est_mnar$simulation_results_sumary )
ggsave("images/simulation_study_results_MNAR_spm_est_100.pdf",
       q_spm_est_mnar,
       width = 20,
       height = 15,
       units = "cm")

format(spm_est_mnar$simulation_results_sumary, digits = 2)
write.csv(format(spm_est_mnar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_spm_param_est_mnar", row.names = FALSE)

### MAR spm estimates ----

spm_est_coverage_spm_mar<-  as.data.table(read.csv(file =  glue::glue("{dir}/SPM_parameters_MAR/100_coverage_spm")))
spm_est_coverage_naive_mar<-  as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MAR/100_coverage_naive")))

spm_est_parameter_est_spm_mar <- as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MAR/100_parameter_est_spm")))
spm_est_parameter_est_naive_mar <-  as.data.table(read.csv(file =glue::glue("{dir}/SPM_parameters_MAR/100_parameter_est_naive")))

spm_est_true_values_mar <- as.data.table(read.csv(file = glue::glue("{dir}/SPM_parameters_MAR/true_values")))
spm_est_mar <- create_sim_data(spm_est_parameter_est_spm_mar, 
                               spm_est_coverage_spm_mar, 
                               spm_est_parameter_est_naive_mar, 
                               spm_est_coverage_naive_mar,
                               spm_est_true_values_mar)

spm_est_mar$simulation_results_sumary

q_spm_est_mar<- plot_sim_results(spm_est_mar$simulation_results, spm_est_mar$simulation_results_sumary )

ggsave("images/simulation_study_results_MAR_spm_est_100.pdf",
       q_spm_est_mar,
       width = 20,
       height = 15,
       units = "cm")

format(spm_est_mar$simulation_results_sumary, digits = 2)
write.csv(format(spm_est_mar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_spm_param_est_mar", row.names = FALSE)


diff_spm_mnar_mar <- abs(spm_est_mnar$simulation_results_sumary$bias_spm) -abs(spm_est_mar$simulation_results_sumary$bias_spm)
diff_naive_mnar_mar <-abs(spm_est_mnar$simulation_results_sumary$bias_naive) -abs(spm_est_mar$simulation_results_sumary$bias_naive)

diff_spm_mnar_mar - diff_naive_mnar_mar
######################################

##################################

### MNAR naive estimates ----

naive_est_coverage_spm_mnar<-  as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MNAR/100_coverage_spm")))
naive_est_coverage_naive_mnar<-  as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MNAR/100_coverage_naive")))

naive_est_parameter_est_spm_mnar <- as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MNAR/100_parameter_est_spm")))
naive_est_parameter_est_naive_mnar <-  as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MNAR/100_parameter_est_naive")))

naive_est_true_values_mnar <- as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MNAR/true_values")))

naive_est_mnar <- create_sim_data(naive_est_parameter_est_spm_mnar, 
                                  naive_est_coverage_spm_mnar, 
                                  naive_est_parameter_est_naive_mnar, 
                                  naive_est_coverage_naive_mnar,
                                  naive_est_true_values_mnar)

naive_est_mnar$simulation_results_sumary

q_naive_est_mnar<- plot_sim_results(naive_est_mnar$simulation_results, naive_est_mnar$simulation_results_sumary )

ggsave("images/simulation_study_results_MNAR_naive_est_100.pdf",
       q_naive_est_mnar,
       width = 20,
       height = 15,
       units = "cm")

format(naive_est_mnar$simulation_results_sumary, digits = 2)
write.csv(format(naive_est_mnar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_naive_est_mnar", row.names = FALSE)

### MAR naive estimates ----

naive_est_coverage_spm_mar<-  as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MAR/100_coverage_spm")))
naive_est_coverage_naive_mar<-  as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MAR/100_coverage_naive")))

naive_est_parameter_est_spm_mar <- as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MAR/100_parameter_est_spm")))
naive_est_parameter_est_naive_mar <-  as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MAR/100_parameter_est_naive")))

naive_est_true_values_mar <- as.data.table(read.csv(file = glue::glue("{dir}/naive_parameters_MAR/true_values")))

naive_est_mar <- create_sim_data(naive_est_parameter_est_spm_mar, 
                                 naive_est_coverage_spm_mar, 
                                 naive_est_parameter_est_naive_mar, 
                                 naive_est_coverage_naive_mar,
                                 naive_est_true_values_mar)

naive_est_mar$simulation_results_sumary

q_naive_est_mar<- plot_sim_results(naive_est_mar$simulation_results, naive_est_mar$simulation_results_sumary )


ggsave("images/simulation_study_results_MAR_naive_est_100.pdf",
       q_naive_est_mar,
       width = 20,
       height = 15,
       units = "cm")

format(naive_est_mar$simulation_results_sumary, digits = 2)
write.csv(format(naive_est_mar$simulation_results_sumary, digits =1 ), file = "results/simulations/summary_naive_est_mar", row.names = FALSE)




