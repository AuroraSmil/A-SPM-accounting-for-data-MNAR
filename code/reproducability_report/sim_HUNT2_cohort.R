
# ### create function ----

create_sim_data <- function(n= 1000){
  set.seed(12345)
  #### Read SPM model and original data #####
  # n = 15000
  # spm_summary_fixed <- as.data.table(read_csv("../../results/SPM/summary_fixed"))
  # setnames(spm_summary_fixed, "X1", "name")
  # 
  # spm_summary_hyper <- as.data.table(read.csv("../../results/SPM/summary_hyper"))
  # 
  # age_table <- as.data.table(read_csv("../../results/SPM/age_effect"))
  # coef_bmi <- as.data.table(read_csv("../../results/reproducability_report/coef_bmi"))
  # coef_bp_2 <- as.data.table(read_csv("../../results/reproducability_report/coef_bp_2"))
  
  spm_summary_fixed <- as.data.table(read_csv("results/SPM/summary_fixed"))
  setnames(spm_summary_fixed, "X1", "name")

  spm_summary_hyper <- as.data.table(read.csv("results/SPM/summary_hyper"))

  age_table <- as.data.table(read_csv("results/SPM/age_effect"))
  coef_bmi <- as.data.table(read_csv("results/reproducability_report/coef_bmi.csv"))
  coef_bp_2 <- as.data.table(read_csv("results/reproducability_report/coef_bp_2.csv"))
  
  
  #n = nrow(data_HUNT23_original)
  
  fake_data = data.table(id = 1:n)
  
  # SEX ----
  p_sex = 1-0.4697
  fake_data[, sex:= rbinom(n, 1, p_sex)]
  
  # AGE ----
  
  mean_age_u <- 74
  std_age_u <- 7
  
  
  mean_age_l <- 44
  std_age_l <- 15
  p = 0.8
  
  split <- rbinom(n, 1, p)
  age_sim_vec <- c(rnorm(n= sum(split), mean = mean_age_l, sd = std_age_l), rnorm(n= n-sum(split), mean = mean_age_u, sd= std_age_u))
  min_age <- min(age_sim_vec)
  while (min_age < 18){
    index <- which(age_sim_vec <18)
    n_new <- length(index)
    #split = rbinom(n_new, 1, p)
    #age_sim_vec_new = c(rnorm(n= sum(split), mean = mean_age_l, sd = std_age_l), rnorm(n= n_new-sum(split), mean = mean_age_u, sd= std_age_u))
    age_sim_vec_new<- rnorm(n= n_new, mean = mean_age_l, sd = std_age_l)
    age_sim_vec[index]<- age_sim_vec_new
    min_age <- min(age_sim_vec)
  }
  
  
  fake_data[, age_2_original_scale := age_sim_vec]
  fake_data[, age_2_original_scale := round(age_sim_vec, 1)]
  fake_data[, age_2 := scale(age_2_original_scale)]
  
  ## BMI ----

  # mean_bmi <-  0 #-0.3
  # std_bmi <-  1
  # fake_data[, bmi_2 := rsnorm(.N, mean = mean_bmi, sd = std_bmi, xi = 1.5)]
  
  mean_bmi <-  0
  std_bmi <-  0.9#2
  p_bmi = 0.5
  #fake_data[, bmi_2 := p_bmi*rsnorm(.N, mean = mean_bmi, sd = std_bmi, xi = 1.5) + (1-p_bmi)*(coef_bmi[1,1][[1]]+ coef_bmi[2,1][[1]] *age_2)]
  fake_data[, bmi_2 := rsnorm(.N, mean = coef_bmi[1,1][[1]]+ coef_bmi[2,1][[1]] *age_2 + coef_bmi[3,1][[1]] *sex, sd = std_bmi, xi = 1.5)]
  
  # BP2 ----
  #fake_data[, bp_2_corr:= rsnorm(.N, 0, 1, 1.5)]# rsnorm(n, 0, 15, 1.5)]
  
  fake_data[, bp_2_eps := rsnorm(.N, coef_bp_2[1,1][[1]] + coef_bp_2[2,1][[1]]*sex + coef_bp_2[3,1][[1]]*age_2 + coef_bp_2[4,1][[1]]*bmi_2-0.15, 0.6, 3)]
  #fake_data[, bp_2_eps := rsnorm(.N, 0, 0.85, 3)]
  #fake_data[, bp_2_noise := rnorm(.N, 0, 0.5)]
  #fake_data[, bp_2_corr:= coef_bp_2[1,1][[1]] + coef_bp_2[2,1][[1]]*sex + coef_bp_2[3,1][[1]]*age_2 + coef_bp_2[4,1][[1]]*bmi_2 +bp_2_eps]# rsnorm(n, 0, 15, 1.5)]
  fake_data[, bp_2_corr:= bp_2_eps]# rsnorm(n, 0, 15, 1.5)]
  
  # BP3 ----
  sigma_eps <- spm_summary_hyper[2,2][[1]]
  fake_data[, eps := rnorm(.N, 0, sigma_eps)]
  fake_data[, bp_3_corr:= spm_summary_fixed[1,2][[1]] + 
              spm_summary_fixed[2,2][[1]]*bp_2_corr + 
              spm_summary_fixed[3,2][[1]]*age_2 + 
              spm_summary_fixed[4,2][[1]]*sex + 
              spm_summary_fixed[5,2][[1]]*bmi_2 +
              eps]
  
  # missing ----

  
  fake_data[, age_round:= round(age_2, 2)]
  
  age_func <- gam(mean ~s(age_round),data=age_table)
  
  fake_data[, age_effekt:= predict.gam(age_func, data.frame(age_round = age_round))]
  c <- spm_summary_hyper[3,2][[1]]
  fake_data[, eta_p:= spm_summary_fixed[6,2][[1]] + 
              spm_summary_fixed[8,2][[1]]*bp_2_corr+ 
              spm_summary_fixed[7,2][[1]]*sex + age_effekt + 
              spm_summary_fixed[9,2][[1]]*bmi_2 + 
              c*eps]
  #fake_data[, eta_p:= coef_spm[6,1] + coef_spm[8,1]*bp_2_corr+ coef_spm[7,1]*sex + -2*cos(age_2) + coef_spm[9,1]*bmi_2 + 0.7*eps]
  
  
  fake_data[, p_i := exp(eta_p)/(exp(eta_p)+1)]
  fake_data[, missing := rbinom(.N, 1, p_i)]
  fake_data[missing==1, bp_3_corr := NA ]
  return(fake_data)
}


create_sim_data(10)
