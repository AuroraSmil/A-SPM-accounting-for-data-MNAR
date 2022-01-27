# Data exploring

library(data.table)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(data.table)
library(latex2exp)
library(grid)
library(gridtext)
library(tikzDevice)

#Sourse basic explore functions
source("code/explore/explore_functions.R") 

# Read data ----
data_HUNT23 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")
data_HUNT34 <- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")

# plot lower and upper bound y axis
ylim_l = 0
ylim_u = 0.7

# Compare missing and pressent participants---- 

q_bp_2 <- smooth_density_compare_missing(data_HUNT23, "bp_2_corr", "$BP_2$", ylim_l, ylim_u, -2.5, 5)
q_age_2 <- smooth_density_compare_missing(data_HUNT23, "age_2", "$Age_2$", ylim_l, ylim_u, -2, 3)
q_bmi_2 <- smooth_density_compare_missing(data_HUNT23, "bmi_2", "$BMI_2$", ylim_l, ylim_u, -4,6)

q_bp_3 <- smooth_density_compare_missing(data_HUNT34, "bp_3_corr", "$BP_3$", ylim_l, ylim_u, -2.5, 5)
q_age_3 <- smooth_density_compare_missing(data_HUNT34, "age_3", "$Age_3$", ylim_l, ylim_u, -2, 3)
q_bmi_3 <- smooth_density_compare_missing(data_HUNT34, "bmi_3", "$BMI_3$", ylim_l, ylim_u, -4,6)

mylegend <- g_legend(q_bp_2)
yleft <- textGrob(expression(paste("Density")), 
                  rot = 90, gp = gpar(fontsize = 15))
density_comp_missing <- grid.arrange(arrangeGrob(
                                                          q_bp_2 + theme(legend.position="none"),
                                                          q_age_2 + theme(legend.position="none"),
                                                          q_bmi_2 + theme(legend.position="none"),
                                                          q_bp_3 + theme(legend.position="none"),
                                                          q_age_3 + theme(legend.position="none"),
                                                          q_bmi_3 + theme(legend.position="none"),
                                                          bottom=mylegend,
                                                          left = yleft,
                                                          nrow=2, 
                                                          ncol = 3),
                                              heights=c(15, 1))


ggsave("images/density_comp_missing.pdf",
       density_comp_missing,
       width = 16,
       height = 15,
       units = "cm")

density_comp_missing_jasa <- grid.arrange(arrangeGrob(
  q_bp_2 + theme(legend.position="none"),
  q_age_2 + theme(legend.position="none"),
  q_bmi_2 + theme(legend.position="none"),
  q_bp_3 + theme(legend.position="none"),
  q_age_3 + theme(legend.position="none"),
  q_bmi_3 + theme(legend.position="none"),
  bottom=mylegend,
  left = yleft,
  nrow=2, 
  ncol = 3),
  heights=c(15, 1))
ggsave("images/density_comp_missing_jasa.pdf",
       density_comp_missing_jasa,
       width = 16,
       height = 7.5,
       units = "cm")


density_comp_missing_jasa_HUNT2 <- grid.arrange(arrangeGrob(
  q_bp_2 + theme(legend.position="none"),
  q_age_2 + theme(legend.position="none"),
  q_bmi_2 + theme(legend.position="none"),
  bottom=mylegend,
  left = yleft,
  nrow=1, 
  ncol = 3),
  heights=c(15, 1))
ggsave("images/density_comp_missing_jasa_HUNT2.pdf",
       density_comp_missing_jasa_HUNT2,
       width = 16,
       height = 5,
       units = "cm")

density_comp_missing_jasa_HUNT3 <- grid.arrange(arrangeGrob(
  q_bp_3 + theme(legend.position="none"),
  q_age_3 + theme(legend.position="none"),
  q_bmi_3 + theme(legend.position="none"),
  bottom=mylegend,
  left = yleft,
  nrow=1, 
  ncol = 3),
  heights=c(15, 1))
ggsave("images/density_comp_missing_jasa_HUNT3.pdf",
       density_comp_missing_jasa_HUNT3,
       width = 16,
       height = 5,
       units = "cm")

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

mean_std_error_bp_3

fit_missing_HUNT3_summary

#Explanatory model HUNT34 ----
formula_bp_4 <- "bp_4_corr~ sex + age_3 + bp_3_corr + bmi_3"

data_relevant_HUNT34 <- as.data.table(data_HUNT34)[, .(sex, age_3, bp_3_corr, bmi_3, bp_4_corr)]
nrow(data_relevant_HUNT34)
data_complete_case_HUNT34 <- na.omit(data_relevant_HUNT34)
nrow(data_complete_case_HUNT34)
fit_lin_bp_4 <- lm(formula_bp_4, data_complete_case_HUNT34)
summary(fit_lin_bp_4)
mean = summary(fit_lin_bp_4)$coefficients

mean_std_error_bp_4 <- data.table(
  variable = c("intercept", "sex", "age_3", "bp_3", "bmi_3"),
  mean = round(fit_lin_bp_4$coefficients, 2),
  std_err = round(coef(summary(fit_lin_bp_4))[, "Std. Error"], 4),
  "Pr(>|t|)" = (coef(summary(fit_lin_bp_4))[, "Pr(>|t|)"])
)

formula_missing_HUNT4 <- "missing~bp_3_corr + sex+ age_3 + bmi_3 "

fit_missing_HUNT4 <- glm(formula_missing_HUNT4, family = "binomial", data = data_HUNT34)

summary(fit_missing_HUNT3)
summary(fit_missing_HUNT4)
fit_missing_HUNT4_summary <- data.table(
  variable = c("intercept", "bp_3", "sex", "age_3", "bmi_3"),
  mean = round(fit_missing_HUNT4$coefficients, 2),
  std_err = round(coef(summary(fit_missing_HUNT4))[, "Std. Error"], 4),
  "Pr(>|z|)" = (coef(summary(fit_missing_HUNT4))[, "Pr(>|z|)"])
)
mean_std_error_bp_4
fit_missing_HUNT4_summary
mean_std_error_bp_3
fit_missing_HUNT3_summary
