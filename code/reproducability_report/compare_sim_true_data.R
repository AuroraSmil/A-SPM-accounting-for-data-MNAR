
library(data.table)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(data.table)
library(latex2exp)
library(tikzDevice)
library(Hmisc)
library(sn)
library(fGarch)
library(INLA)
library(brinla)
library(tidyverse)
library(data.table)
library(latex2exp)
library(lubridate)
library(glue)
library(readr)
library(corrplot)
library(grid)
library(gridtext)
library(mgcv)

source("code/explore/explore_functions.R")
source("code/simulation_study/model_functions.R")
source("code/reproducability_report/sim_HUNT2_cohort.R", local = knitr::knit_global())
# Compare simulated data to true data
# 
n = 64385#1000
fake_data <- create_sim_data(n)
write_csv(fake_data, glue("code/reproducability_report/sim_data_{n}.csv"))
data_HUNT23_scaled <- as.data.table(readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")[, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing)])
data_HUNT23 <- as.data.table(readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_adjusted_meds.RData")[, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing)])

# Plot an display fake data ----
setnames(fake_data, c("age_2", "bp_2_corr", "bmi_2", "bp_3_corr", "missing"), 
         c("age", "BP_2", "BMI", "BP_3", "m"))
setnames(data_HUNT23_scaled, c("age_2", "bp_2_corr", "bmi_2", "bp_3_corr", "missing"), 
         c("age", "BP_2", "BMI", "BP_3", "m"))


corr_true<- rcorr(as.matrix(fake_data[, .( sex, age, BP_2, BMI, BP_3, m)]))[1]$r
corr_sim <- rcorr(as.matrix(data_HUNT23_scaled[, .( sex, age, BP_2, BMI, BP_3, m)]))[1]$r


corrplot_true <- { # Prepare the Corrplot 
  corrplot(corr_true,  
           type = "lower", 
           method = "circle",
           tl.cex=1,
           #title = "Regional Factor Correlation Matrix over history", 
           mar = c(0,0,1,0), 
           number.cex = 0.5, 
           number.digits = 2);
  # Call the recordPlot() function to record the plot
  recordPlot()
}
corrplot_true

write_csv(as.data.frame(corr_true), "results/reproducability_report/correlation_true.csv")

# In case if you want to save the image using ggsave
# replayPlot basically prints the plot.
ggsave(filename = "images/corrplot_true.pdf", 
       plot = replayPlot(corrplot_true))



corrplot_sim <- { # Prepare the Corrplot 
  corrplot(corr_sim,  
           type = "lower", 
           method = "circle",
           #title = "Regional Factor Correlation Matrix over history", 
           mar = c(0,0,1,0), 
           number.cex = 0.5, 
           number.digits = 2);
  # Call the recordPlot() function to record the plot
  recordPlot()
}
write_csv(as.data.frame(corr_sim), "results/reproducability_report/correlation_sim.csv")
# In case if you want to save the image using ggsave
# replayPlot basically prints the plot.
ggsave(filename = "images/corrplot_sim.pdf", 
       plot = replayPlot(corrplot_sim))



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

q <- ggplot(data = data_HUNT23_scaled, aes(x = age))
q <- q + geom_density(aes(color = "Observed data", linetype = "Observed data"))
q <- q + geom_density(data = fake_data, aes(x = age, color = "Simulated data", linetype = "Simulated data"))
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + labs(color  = "Guide name", linetype = "Guide name")
q <- q + xlab(unname(unname(TeX(c("Age"))))) + ylab(element_blank())
q_age <- q


q <- ggplot(data = data_HUNT23_scaled, aes(x = BMI))
q <- q + geom_density(aes(color = "Observed data", linetype = "Observed data"))
q <- q + geom_density(data = fake_data, aes(x = BMI, color = "Simulated data", linetype = "Simulated data"))
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + labs(color  = "Guide name", linetype = "Guide name")
q <- q + xlab(unname(unname(TeX(c("BMI"))))) + ylab(element_blank())
q_bmi <- q



q <- ggplot(data = data_HUNT23_scaled, aes(x = BP_2))
q <- q + geom_density(aes(color = "Observed data", linetype = "Observed data"))
q <- q + geom_density(data = fake_data, aes(x = BP_2, color = "Simulated data", linetype = "Simulated data"))
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + labs(color  = "Guide name", linetype = "Guide name")
q <- q + xlab(unname(unname(TeX(c("BP_2"))))) + ylab(element_blank())
q_bp_2 <- q



q <- ggplot(data = data_HUNT23_scaled, aes(x = BP_3))
q <- q + geom_density(aes(color = "Observed data", linetype = "Observed data"))
q <- q + geom_density(data = fake_data, aes(x = BP_3, color = "Simulated data", linetype = "Simulated data"))
q <- q + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + labs(color  = "Guide name", linetype = "Guide name")
q <- q + xlab(unname(unname(TeX(c("BP_3"))))) + ylab(element_blank())
q_bp_3<- q



mylegend <- g_legend(q_age)
yleft <- textGrob(expression(paste("Density")),
                  rot = 90, gp = gpar(fontsize = 15))
hist_comp_sim <- grid.arrange(arrangeGrob(
  q_age + theme(legend.position="none"),
  q_bmi + theme(legend.position="none"),
  q_bp_2 + theme(legend.position="none"),
  q_bp_3 + theme(legend.position="none"),
  bottom=mylegend,
  left = yleft,
  nrow=2,
  ncol = 2),
  heights=c(15, 1))

ggsave("images/hist_comp_true_sim.pdf",
       hist_comp_sim,
       width = 15,
       height = 15,
       units = "cm")

ylim_l = 0
ylim_u = 0.6


# q_bp_2 <- smooth_density_compare_missing(data_HUNT23_scaled, "bp_2_corr", "$BP_2$", ylim_l, ylim_u)
# q_age_2 <- smooth_density_compare_missing(data_HUNT23_scaled, "age_2", "$Age_2$", ylim_l, ylim_u)
# q_bmi_2 <- smooth_density_compare_missing(data_HUNT23_scaled, "bmi_2", "$BMI_2$", ylim_l, ylim_u)
# 
# mylegend <- g_legend(q_bp_2)
# yleft <- textGrob(expression(paste("Density")), 
#                   rot = 90, gp = gpar(fontsize = 15))
# hist_comp_missing <- grid.arrange(arrangeGrob(
#   q_bp_2 + theme(legend.position="none"),
#   q_age_2 + theme(legend.position="none"),
#   q_bmi_2 + theme(legend.position="none"),
#   bottom=mylegend,
#   left = yleft,
#   nrow=1, 
#   ncol = 3),
#   heights=c(15, 1))
# 
# ggsave("images/hist_comp_HUNT23_missing.pdf",
#        hist_comp_missing,
#        width = 15,
#        height = 15,
#        units = "cm")

