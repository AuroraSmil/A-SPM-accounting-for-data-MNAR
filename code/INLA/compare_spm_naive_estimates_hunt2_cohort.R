# library(cowplot)
# library(gridExtra)
# library(grid)
# library(gridtext)
# library(INLA)
# library(data.table)
# library(ggplot2)
# library(tikzDevice)
# library(latex2exp)
# 
# library(INLA)
# library(brinla)
# #INLA:::inla.dynload.workaround()
# library(tidyverse)
# library(data.table)
# library(latex2exp)
# library(tikzDevice)
# library(cowplot)
# library(gridExtra)
# library(grid)
# library(gridtext)
# source("code/simulation_study/model_functions.R")
# source("code/explore/explore_functions.R")


compare_spm_naive_estimates <- function(dir, cohort_number = 2){
  
  # read results ----
  
  if(cohort_number == 2){
    inla_spm <- readRDS(file = glue::glue("{dir}spm_inla_HUNT23.RData"))
    inla_naive_bp <- readRDS(file = glue::glue("{dir}inla_naive_bp.RData" ))
    inla_naive_m <- readRDS(file =  glue::glue("{dir}inla_naive_m.RData"))
    
  }else if(cohort_number == 3){
    inla_spm <- readRDS(file = glue::glue("{dir}spm_inla_HUNT34.RData"))
    inla_naive_bp <- readRDS(file = glue::glue("{dir}inla_naive_bp_HUNT34.RData" ))
    inla_naive_m <- readRDS(file =  glue::glue("{dir}inla_naive_m_HUNT34.RData"))
    
  }else{
    print("unknown cohort")
    break
  }

  bri.hyperpar.summary(inla_spm)
  bri.hyperpar.summary(inla_naive_bp)
  bri.hyperpar.summary(inla_naive_m)
  
  #Constructt a vector to store plots
  
  plots_naive_spm <- vector("list", 
                            length = (length(names(inla_naive_bp$marginals.fixed)) 
                                      +length(names(inla_naive_m$marginals.fixed))))
  
  # Vector with names for the plots
  naive_fixed_names <- c(names(inla_naive_bp$marginals.fixed), names(inla_naive_m$marginals.fixed))
  naive_fixed_names_tex <- c("$\\alpha_0$", "$\\alpha_{BP}$", "$\\alpha_{age}$" ,  "$\\alpha_{sex}$"  , "$\\alpha_{BMI}$", 
                             "$\\beta_0$", "$\\beta_{BP}$",  "$\\beta_{sex}$"  , "$\\beta_{BMI}$" )
  
  #Plots for the fixed effects related to bp -----
  for(i in 1:length(inla_naive_bp$marginals.fixed)){
    data_plot_n <- as.data.table(inla_naive_bp$marginals.fixed[i])
    data_plot_spm <- as.data.table(inla_spm$marginals.fixed[i])
    setnames(data_plot_n, c("x", "y_n"))
    setnames(data_plot_spm, c("x", "y_spm"))
    q <- ggplot(data = data_plot_n, aes(x = x))
    q <- q + geom_line(aes(y = y_n, colour = "Naive", linetype = "Naive"))
    q <- q + geom_line(data = data_plot_spm, aes(y = y_spm, colour = "SPM", linetype = "SPM"))
    q <- q + xlab(unname(TeX(c(glue::glue("{naive_fixed_names_tex[i]}")))))
    q <- q + ylab("Density")
    q <- q + ylab(element_blank())
    q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
    q <- q + labs(colour  = "Guide name", linetype = "Guide name")
    q <- q + scale_linetype_manual(values = c("solid","dashed"))
    q
    name <- naive_fixed_names[[i]]
    plots_naive_spm[[name]] <- q
  }
  
  #Plots for the fixed effects related to the drop out prosess ----
  for(i in 1:length(inla_naive_m$marginals.fixed)){
    data_plot_n <- as.data.table(inla_naive_m$marginals.fixed[i])
    data_plot_spm <- as.data.table(inla_spm$marginals.fixed[i+ length(inla_naive_bp$marginals.fixed)])
    setnames(data_plot_n, c("x", "y_n"))
    setnames(data_plot_spm, c("x", "y_spm"))
    q <- ggplot(data = data_plot_n, aes(x = x))
    q <- q + geom_line(aes(y = y_n, colour = "Naive", linetype = "Naive"))
    q <- q + geom_line(data = data_plot_spm, aes(y = y_spm, colour = "SPM", linetype = "SPM"))
    q <- q + xlab(unname(TeX(c(glue::glue("{naive_fixed_names_tex[i + length(inla_naive_bp$marginals.fixed)]}")))))
    q <- q + ylab("Density")
    q <- q + ylab(element_blank())
    q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
    q <- q + labs(colour  = "Guide name", linetype = "Guide name")
    q <- q + scale_linetype_manual(values = c("solid","dashed"))
    
    name <- naive_fixed_names[[i+length(inla_naive_bp$marginals.fixed)]]
    plots_naive_spm[[name]] <- q
  }
  
  # Precision data for the hyperparameters and convert it to standard deviation
  precision_eps_n_bp<- inla_naive_bp$marginals.hyperpar$`Precision for y_eps1`
  sigma_eps_n_bp<- inla.tmarginal(function(x) 1/sqrt((x)), precision_eps_n_bp)
  precision_eps_spm<- inla_spm$marginals.hyperpar$`Precision for y_eps1`
  sigma_eps_spm <- inla.tmarginal(function(x) 1/sqrt((x)), precision_eps_spm)
  
  precision_eps_n_age<- inla_naive_m$marginals.hyperpar$`Precision for m_AGE`
  sigma_eps_n_age<- inla.tmarginal(function(x) 1/sqrt(x), precision_eps_n_age)
  precision_eps_spm_age<- inla_spm$marginals.hyperpar$`Precision for m_AGE`
  sigma_eps_spm_age <- inla.tmarginal(function(x) 1/sqrt(x), precision_eps_spm_age)
  
  
  # Plot sigma_epsilon ----
  q <- ggplot(data = as.data.frame(sigma_eps_n_bp), aes(x = x, y = y))
  q <- q + geom_line(aes( colour = "Naive",linetype = "Naive"))
  q <- q + geom_line(data = as.data.frame(sigma_eps_spm), aes(y = y, colour = "SPM", linetype = "SPM"))
  q <- q + xlab(unname(TeX(c("$\\sigma_{\\epsilon}"))))
  q <- q + ylab("Density")
  q <- q + ylab(element_blank())
  q <- q + theme(legend.title = element_blank())
  q <- q + labs(colour  = "Guide name", linetype = "Guide name")
  q <- q + scale_linetype_manual(values = c("solid","dashed"))
  
  
  plots_naive_spm$sigma_eps <- q
  
  # Plot sigma_age ----
  q <- ggplot(data = as.data.frame(sigma_eps_n_age), aes(x = x, y = y))
  q <- q + geom_line(aes( colour = "Naive", linetype = "Naive"))
  q <- q + geom_line(data = as.data.frame(sigma_eps_spm_age), aes(y = y, colour = "SPM", linetype = "SPM"))
  q <- q + xlab(unname(TeX(c("$\\sigma_{age}"))))
  q <- q + ylab("Density") 
  q <- q + ylab(element_blank())
  q <- q + theme(legend.title = element_blank())
  q <- q + labs(colour  = "Guide name", linetype = "Guide name")
  q <- q + scale_linetype_manual(values = c("solid","dashed"))
  
  
  plots_naive_spm$sigma_age <- q
  
  # Plot assosiation parameter c ----
  c<- inla_spm$marginals.hyperpar$`Beta for m_eps1`
  
  q <- ggplot(data = as.data.frame(c), aes(x = x))
  q <- q + geom_line(aes(y = y, colour = "SPM"), linetype = "dashed")
  q <- q + xlab(unname(TeX(c("c"))))
  q <- q + ylab("Density")
  q <- q + ylab(element_blank()) 
  q <- q + theme(legend.title = element_blank())
  q <- q + scale_color_manual(values=c("#00BFC4"))

  
  plots_naive_spm$c<- q
  
  # Plot the  ageeffect ----
  
  data_random_age_inla <-  inla_spm$summary.random$m_AGE
  data_random_age_naive <-  inla_naive_m$summary.random$m_AGE
  q <- ggplot(data = data_random_age_inla, aes(x = ID))
  q <- q + geom_line(aes(y = mean, colour = "SPM", linetype = "SPM"))
  q <- q + geom_ribbon(aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "SPM"), alpha = 0.5)
  q <- q + geom_line(data = data_random_age_naive, aes(y = mean, colour = "Naive", linetype = "Naive"))
  q <- q + geom_ribbon(data = data_random_age_naive, aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "naive"), alpha = 0.5)
  q <- q + labs(colour  = "Guide name", linetype = "Guide name")
  q <- q + scale_linetype_manual(values = c("solid","dashed"))
  q_age <- q + xlab(TeX("Age")) + ylab("Age effect") +guides(fill = FALSE) + theme( legend.position = "bottom", legend.title = element_blank())

  plots_naive_spm$age_effect<- q_age
  
  
  # Display all graphs in one figure
  
  mylegend<-g_legend(plots_naive_spm$alpha_0)
  yleft <- textGrob(expression(paste("Density")), 
                    rot = 90, gp = gpar(fontsize = 15))
  
  post_dens_coef <- grid.arrange(arrangeGrob(plots_naive_spm$alpha_0 + theme(legend.position="none"),
                                             plots_naive_spm$y_BP1 + theme(legend.position="none"),
                                             plots_naive_spm$y_BMI + theme(legend.position="none"),
                                             plots_naive_spm$y_AGE + theme(legend.position="none"),
                                             plots_naive_spm$y_SEX + theme(legend.position="none"),
                                             plots_naive_spm$sigma_eps + theme(legend.position="none"),
                                             plots_naive_spm$c + theme(legend.position="none"),
                                             plots_naive_spm$beta_0 + theme(legend.position="none"),
                                             plots_naive_spm$m_BP1 + theme(legend.position="none"),
                                             plots_naive_spm$m_BMI + theme(legend.position="none"),
                                             plots_naive_spm$m_SEX + theme(legend.position="none"),
                                             plots_naive_spm$sigma_age + theme(legend.position="none"), 
                                             ncol= 2),
                                 nrow = 2, 
                                 left = yleft,
                                 mylegend,
                                 heights=c(12, 1))
  ggsave(glue::glue("images/post_dens_coef_scaled_hunt{cohort_number}_param_est.pdf"),
         post_dens_coef,
         width = 15,
         height = 20,
         units = "cm")
  ggsave(glue::glue("images/age_effect_hunt{cohort_number}.pdf"),
         q_age,
         width = 15,
         height = 7.5,
         units = "cm")
  
  # Store summary results ----
  
  if(cohort_number == 2){
    inla_spm <- readRDS( file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT23.RData")
    summary_fixed <- as.data.frame(round(inla_spm$summary.fixed,3))
    summary_hyper <- as.data.frame(round(bri.hyperpar.summary(inla_spm),3))
    age_table <- as.data.table(inla_spm$summary.random$m_AGE[, c("ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")])
    age_table[, age_round := round(ID, 2)]
    
    
    write.csv(summary_fixed, "results/SPM/summary_fixed")
    write.csv(summary_hyper, "results/SPM/summary_hyper")
    write.csv(age_table, "results/SPM/age_effect", row.names = FALSE)
    
  }else if (cohort_number == 3){
    inla_spm <- readRDS( file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT34.RData")
    summary_fixed <- as.data.frame(round(inla_spm$summary.fixed,3))
    summary_hyper <- as.data.frame(round(bri.hyperpar.summary(inla_spm),3))
    age_table <- as.data.table(inla_spm$summary.random$m_AGE[, c("ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")])
    age_table[, age_round := round(ID, 2)]
    
    
    write.csv(summary_fixed, "results/SPM/summary_fixed_HUNT34")
    write.csv(summary_hyper, "results/SPM/summary_hyper_HUNT34")
    write.csv(age_table, "results/SPM/age_effect_HUNT34", row.names = FALSE)
    
  }
  

}


