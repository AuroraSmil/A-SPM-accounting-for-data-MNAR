# Data exploring

library(data.table)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(data.table)
library(latex2exp)
library(tikzDevice)

#Sourse basic explore functions
source("code/explore/explore_functions.R") 

data_HUNT12 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT12_scaled.RData")
#data_HUNT12 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT12_adjusted_meds.RData")

data_HUNT23 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled_HUNT12.RData")
#data_HUNT23 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_adjusted_meds.RData")

ylim_l = 0
ylim_u = 0.75



# Look for correlation between missingness and non missing---- 
q_bp_1 <- smooth_density_compare_missing(data_HUNT12, "bp_1_corr", "$BP_1$", ylim_l, ylim_u)
q_age_1 <- smooth_density_compare_missing(data_HUNT12, "age_1", "$Age_1$", ylim_l, ylim_u)
q_bmi_1 <- smooth_density_compare_missing(data_HUNT12, "bmi_1", "$BMI_1$", ylim_l, ylim_u)

q_bp_2 <- smooth_density_compare_missing(data_HUNT23, "bp_2_corr", "$BP_2$", ylim_l, ylim_u)
q_age_2 <- smooth_density_compare_missing(data_HUNT23, "age_2", "$Age_2$", ylim_l, ylim_u)
q_bmi_2 <- smooth_density_compare_missing(data_HUNT23, "bmi_2", "$BMI_2$", ylim_l, ylim_u)

mylegend <- g_legend(q_bp_1)
yleft <- textGrob(expression(paste("Density")), 
                  rot = 90, gp = gpar(fontsize = 15))
hist_comp_missing_HUNT1_HUNT2 <- grid.arrange(arrangeGrob(q_bp_1 + theme(legend.position="none"),
                                               
                                               q_age_1 + theme(legend.position="none"),
                                               q_bmi_1 + theme(legend.position="none"),
                                               q_bp_2 + theme(legend.position="none"),
                                               q_age_2 + theme(legend.position="none"),
                                               
                                               q_bmi_2 + theme(legend.position="none"),
                                               
                                               bottom=mylegend,
                                               
                                               nrow=2, 
                                               ncol = 3),
                                   heights=c(15, 1))

hist_comp_missing_HUNT1_HUNT2
ggsave("images/hist_comp_missing_HUNT1_HUNT2.pdf",
       hist_comp_missing_HUNT1_HUNT2,
       width = 15,
       height = 15,
       units = "cm")


#Explanatory model HUNT12
formula_bp_2 <- "bp_2_corr~ sex + age_1 + bp_1_corr + bmi_1"

data_relevant_HUNT12 <- as.data.table(data_HUNT12)[, .(sex, age_1, bp_1_corr, bmi_1, bp_2_corr)]
nrow(data_relevant_HUNT12)
data_complete_case_HUNT12 <- na.omit(data_relevant_HUNT12)
nrow(data_complete_case_HUNT12)
fit_lin_bp_2 <- lm(formula_bp_2, data_complete_case_HUNT12)
summary(fit_lin_bp_2)
mean = summary(fit_lin_bp_2)$coefficients

mean_std_error_bp_2 <- data.table(
  variable = c("intercept", "sex", "age_1", "bp_1", "bmi_1"),
  mean = round(fit_lin_bp_2$coefficients, 2),
  std_err = round(coef(summary(fit_lin_bp_2))[, "Std. Error"], 4),
  "Pr(>|t|)" = (coef(summary(fit_lin_bp_2))[, "Pr(>|t|)"])
)


formula_missing_HUNT2 <- "missing~bp_1_corr + sex+ age_1 + bmi_1 "

fit_missing_HUNT2 <- glm(formula_missing_HUNT2, family = "binomial", data = data_HUNT12)



fit_missing_HUNT2_summary <- data.table(
  variable = c("intercept", "bp_1", "sex", "age_1", "bmi_1"),
  mean = round(fit_missing_HUNT2$coefficients, 2),
  std_err = round(coef(summary(fit_missing_HUNT2))[, "Std. Error"], 4),
  "Pr(>|z|)" = (coef(summary(fit_missing_HUNT2))[, "Pr(>|z|)"])
)


#Explanatory model HUNT23
formula_bp_3 <- "bp_3_corr~ sex + age_2 + bp_2_corr + bmi_2"

data_relevant_HUNT23 <- as.data.table(data_HUNT23)[, .(sex, age_2, bp_2_corr, bmi_2, bp_3_corr)]
nrow(data_relevant_HUNT23)
data_complete_case_HUNT23 <- na.omit(data_relevant_HUNT23)
nrow(data_complete_case_HUNT23)
fit_lin_bp_3 <- lm(formula_bp_3, data_complete_case_HUNT23)
summary(fit_lin_bp_3)
mean = summary(fit_lin_bp_3)$coefficients

mean_std_error_bp_3 <- data.table(
  variable = c("intercept", "sex", "age_2", "bp_2", "bmi_2"),
  mean = round(fit_lin_bp_3$coefficients, 2),
  std_err = round(coef(summary(fit_lin_bp_3))[, "Std. Error"], 4),
  "Pr(>|t|)" = (coef(summary(fit_lin_bp_3))[, "Pr(>|t|)"])
)


formula_missing_HUNT3 <- "missing~bp_2_corr + sex+ age_2 + bmi_2 "

fit_missing_HUNT3 <- glm(formula_missing_HUNT3, family = "binomial", data = data_HUNT23)



fit_missing_HUNT3_summary <- data.table(
  variable = c("intercept", "bp_2", "sex", "age_2", "bmi_2"),
  mean = round(fit_missing_HUNT3$coefficients, 2),
  std_err = round(coef(summary(fit_missing_HUNT3))[, "Std. Error"], 4),
  "Pr(>|z|)" = (coef(summary(fit_missing_HUNT3))[, "Pr(>|z|)"])
)


mean_std_error_bp_2
mean_std_error_bp_3

fit_missing_HUNT2_summary
fit_missing_HUNT3_summary
