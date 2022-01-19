library(INLA)
library(brinla)
#INLA:::inla.dynload.workaround()
library(tidyverse)
library(data.table)
library(latex2exp)
library(tikzDevice)
library(cowplot)
library(gridExtra)
library(grid)
library(gridtext)
library(cowplot)
library(ggplot2)
library(glue)
library(dplyr)
library(DescTools)
library(scoringRules)
library(mgcv)



source("code/simulation_study/model_functions.R")
source("code/explore/explore_functions.R")
source("code/INLA/compare_spm_naive_estimates_hunt2_cohort.R")
source("code/INLA/assosiation_effect.R")
# Set directory for storing the SPM and naive model
dir <- "/home/aurorach/data_HUNT_aurora/master/"
# Read data ---- 
data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")
data_small <- as.data.table(data)#[1:10000])
n <- nrow(data_small)
data <- data_small

sigma2 = 0.001 #Fixed to a neglectable value

# Run SPM ---- 

inla_spm <- run_spm(data, verbose = T, compute_all = TRUE)


saveRDS(inla_spm, file = glue::glue{"{dir}spm_inla_HUNT23.RData"})
summary(inla_spm)

# Run naive -----
inla_naive_bp <- run_naive_bp(data, verbose = T, compute_all = TRUE)
saveRDS(inla_naive_bp, file = glue::glue{"{dir}inla_naive_bp.RData"})

inla_naive_m <- run_naive_m(data, verbose = T, compute_all = TRUE)
saveRDS(inla_naive_m, file = glue::glue{"{dir}inla_naive_m.RData"})

summary(inla_naive_bp)
summary(inla_naive_m)

# Compare SPM and naive model estimates ----
# Can be run without the above code, just remember to load the packages.

compare_spm_naive_estimates(dir)

# Explore assosiation effect ----

assosiatoin_effect(dir)


