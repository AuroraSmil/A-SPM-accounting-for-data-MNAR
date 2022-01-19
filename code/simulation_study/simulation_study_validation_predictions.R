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



simulation_study_validation_predictions_func <- function(model = "SPM", n_sim = 100, n_samples = 300){
  # Create directory for storing results 
  today <- today()
  if(!dir.exists(glue("results/simulations/{today}"))){
    dir.create(glue("results/simulations/{today}"))
  }
  if(!dir.exists(glue("results/simulations/{today}/{model}"))){
    dir.create(glue("results/simulations/{today}/{model}"))
  }
  
  dir <- glue("results/simulations/{today}/{model}")
  
  # Store simulation results 
  crps_vec <- vector(mode = "list", length = n_sim)
  brier_vec <- vector(mode = "list", length = n_sim)
  
  
  if(model == "SPM"){
    inla_spm <-readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_sim_HUNT23.RData")
  }else if (model == "naive"){
    naive_bp <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/inla_naive_bp_sim_HUNT23.RData")
    naive_m <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/inla_naive_m_sim_HUNT23.RData")
  }else{
    stop("Model is wrong")
  }
  
  
  # simulations starts -----
  start_overall <- Sys.time()
  nr_failes <- 0
  i = 1
  print("Starting simulatons")
  while(i <=n_sim){
    print(i)
    #create data sim
    start_time <- Sys.time()
    data_sim <- as.data.table(read.csv2(file = glue::glue("/home/aurorach/data_HUNT_aurora/master/sim_data_validation/data_sim_{i}")))
    
    # Validation function
    
    if(model == "SPM"){
      validation <- try(validation_func_spm(spm = inla_spm,
                                            data = data_sim, 
                                            n_samples = n_samples))
    }else if (model == "naive"){
      validation <- try(validation_func_naive(naive_bp = naive_bp,
                                              naive_m = naive_m,
                                              data = data_sim, 
                                              n_samples = n_samples))
    }else{
      stop("Model is wrong")
    }
    
    if(class(validation) == "try-error"){
      print("validation func crashed")
      nr_failes <- nr_failes + 1
    }else{
      predictions_m <- validation$predictions_m
      predictions_bp <- validation$predictions_bp
      
      # Compute CRPS and Brier
      crps <- crps_func_general(predictions_bp = predictions_bp,  data_val = data_sim, n_samples = n_samples)
      
      crps[, sim_id := i]
      crps_vec[[i]] <- crps
      print(crps)
      brier <- brier_func_general(predictions_m = predictions_m, data_val = data_sim)
      brier[, sim_id := i]
      brier_vec[[i]] <- brier
      print(brier)
    }
    
    end_time <- Sys.time()
    print(glue::glue("Time round {i}: {end_time- start_time}"))
    print(glue::glue("Nr failes so far: {nr_failes}"))
    
    if(!i%%50){
      crps_table <- rbindlist(crps_vec)
      brier_table <- rbindlist(brier_vec)
      write.csv(crps_table, file = glue("{dir}/{i}_crps_sim_model"), row.names = FALSE)
      write.csv(brier_table, file = glue("{dir}/{i}_brier_sim_model"), row.names = FALSE)
    }
    if(!class(validation)=="try-error"){
      i <- i + 1
    }
  }
  
  # Merge all results and store them
  
  end_overall <- Sys.time()
  print(glue::glue("Time round {i}: {end_overall- start_overall}"))
  crps_table <- rbindlist(crps_vec)
  brier_table <- rbindlist(brier_vec)
  write.csv(crps_table, file = glue("{dir}/{n_sim}_crps_sim_model"), row.names = FALSE)
  write.csv(brier_table, file = glue("{dir}/{n_sim}_brier_sim_model"), row.names = FALSE)
  
}


simulation_study_validation_predictions_func(model = "SPM")
simulation_study_validation_predictions_func(model = "naive")
