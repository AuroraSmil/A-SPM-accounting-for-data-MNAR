library(INLA)
#library(inlabru)
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

## Read the models and data ----

inla_spm <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/spm_inla_HUNT23.RData")
inla_naive_bp <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/inla_naive_bp.RData" )
inla_naive_m <- readRDS(file =  "/home/aurorach/data_HUNT_aurora/master/inla_naive_m.RData")

data_val<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")
#data_val <- data_val[1:10, ]

# Run validation function to obtain predictions 
n_samples <- 300
validation <- validation_func(spm = inla_spm, naive_bp = inla_naive_bp, naive_m = inla_naive_m, data = data_val, n_samples = n_samples)
predictions_m_n <- validation$predictions_m_n
predictions_m_spm <- validation$predictions_m_spm
predictions_bp_n <- validation$predictions_bp_n
predictions_bp_spm <- validation$predictions_bp_spm

# SPM ----
# Melt data table 
predictions_bp_long_spm <- melt(predictions_bp_spm, id.vars= "id")
setnames(predictions_bp_long_spm, c("value", "variable"), c("bp_pred", "row_id"))

predictions_m_long_spm <- melt(predictions_m_spm, id.vars= c("id"))
setnames(predictions_m_long_spm, c("value", "variable"), c("m_pred", "row_id"))


# Merge and aggregate data
predictions_long_spm <- as.data.table(merge(predictions_bp_long_spm, predictions_m_long_spm, by = c("row_id", "id")))

predictions_agg_spm <- predictions_long_spm[, .(q500_bp = q500(bp_pred),
                                                mean_bp = mean(bp_pred),
                                         q025_bp = q025(bp_pred), 
                                         q975_bp = q975(bp_pred),
                                         q500_m = q500(m_pred),
                                         mean_m = mean(m_pred),
                                         q025_m = q025(m_pred),
                                         q975_m = q975(m_pred)
                                        ), 
                                    keyby = .(id)]
n_val <- nrow(data_val)
data_val[, id:= 1:n_val]
predictions_agg_spm[data_val, on = "id", bp_corr:= bp_4_corr]
predictions_agg_spm[data_val, on = "id", missing:= missing]
predictions_agg_spm[, residuals_bp := q500_bp - bp_corr]
predictions_agg_spm[, diff_m := q500_m - missing]

mean(abs(na.omit(predictions_agg_spm$mean_bp- predictions_agg_spm$bp_corr)))

# Naive ---- 
# Melt merge and aggregate
predictions_bp_long_n <- melt(predictions_bp_n, id.vars= "id")
setnames(predictions_bp_long_n, c("value", "variable"), c("bp_pred", "row_id"))

predictions_m_long_n <- melt(predictions_m_n, id.vars= c("id"))
setnames(predictions_m_long_n, c("value", "variable"), c("m_pred", "row_id"))

predictions_long_n <- as.data.table(merge(predictions_bp_long_n, predictions_m_long_n, by = c("row_id", "id")))

predictions_agg_n <- predictions_long_n[, .(q500_bp = q500(bp_pred),
                                            mean_bp = mean(bp_pred),
                                                q025_bp = q025(bp_pred), 
                                                q975_bp = q975(bp_pred),
                                            mean_m = mean(m_pred),
                                                q500_m = q500(m_pred),
                                                q025_m = q025(m_pred),
                                                q975_m = q975(m_pred)
), 
keyby = .(id)]

data_val[, id:= 1:n_val]
predictions_agg_n[data_val, on = "id", bp_corr:= bp_4_corr]
predictions_agg_n[data_val, on = "id", missing:= missing]
predictions_agg_n[, residuals_bp := bp_corr - q500_bp]
predictions_agg_n[, diff_m := ( missing- q500_m)]

mean(abs(na.omit(predictions_agg_spm$mean_bp- predictions_agg_spm$bp_corr)))
mean(abs(na.omit(predictions_agg_n$mean_bp - predictions_agg_n$bp_corr)))


###### CRPS ----

crps <- crps_func(predictions_bp_spm = predictions_bp_spm, predictions_bp_n = predictions_bp_n, data_val = data_val, n_samples = n_samples)
crps


###### BRIER -----

brier <- brier_func(predictions_m_spm = predictions_m_spm, predictions_m_n = predictions_m_n, data_val = data_val)
brier

##### compare grouped by missing#####
predictions_agg_spm[, model := "SPM posterior predictive mean"]
predictions_agg_n[, model := "Naive posterior predicitve mean"]
predictions_agg <- rbindlist(list(predictions_agg_spm, predictions_agg_n), use.names = TRUE)
predictions_agg[missing == 0, missing_name := "Present"]
predictions_agg[missing == 1, missing_name := "Missing"]

q <- ggplot(predictions_agg, aes(x = mean_bp, group = model, colour = model, linetype = model))
q <- q + geom_density()
q <- q + geom_density(aes(x = bp_corr, colour = "Observed", linetype = "Observed"))
q <- q + scale_color_manual(values = c("#F8766D", "black","#00BFC4"))
q <- q + xlab(TeX("$BP_4$")) + ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + facet_wrap(vars(missing_name))
q



ggsave("images/validation_pred_bp_smooth.pdf",
       width = 17,
       height = 10,
       units = "cm")

predictions_agg_mean <- predictions_agg[, .(mean_mean_bp = mean(mean_bp), 
                                            mean_mean_m = mean(mean_m)),
                                        keyby  = .(missing, missing_name)]


q <- ggplot(predictions_agg, aes(x = mean_m, group = model, colour = model, linetype = model))
q <- q + geom_density()
#q <- q + xlab("Probability of dropping out") 
q <- q + xlab(TeX("p_4")) 
q <- q + ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
q <- q + scale_color_manual(values = c("#F8766D","#00BFC4", "black"))
q <- q+ facet_wrap(vars(missing_name))
q
ggsave("images/validation_pred_missingness_smooth.pdf",
       width = 17,
       height = 10,
       units = "cm")


##### Construct simulated individuals ----
set.seed(1234)
data_fake <- data.table(
  id = c(1,2,3),
  bp_3_corr = c(-2, 0, 2),
  age_3 = c(-1.5, 0, 1.5),
  bmi_3 = c(-2, 0, 2),
  sex = c(0,0,0),
  name_id = c("Id = 1", "Id = 2", "Id = 3"),
  name = c("Participant 1 ", "Participant 2", "Participant 3")
)

# data_fake <- data.table(
#   id = c(1,2,3,4,5,6),
#   bp_3_corr = c(-2, 0, 2, -2, 0, 2),
#   age_3 = c(-1.5, 0, 1.5, -1.5, 0, 1.5),
#   bmi_3 = c(-2, 0, 2, -2, 0, 2),
#   sex = c(0,0,0, 1, 1,1),
#   name_id = c("ID = 1", "Id = 2", "Id = 3"),
#   name = c("Current systolic bloodpressure = -1.7 ", "Current systolic bloodpressure = 0", "Current systolic bloodpressure = 1.7")
# )

validation_fake <- validation_func(spm = inla_spm, naive_bp = inla_naive_bp, naive_m = inla_naive_m, data_fake, n_samples = n_samples)

fake_predictions_bp_spm <- validation_fake$predictions_bp_spm
fake_predictions_bp_long_spm <- melt(fake_predictions_bp_spm, id.vars= c("id"))
setnames(fake_predictions_bp_long_spm, c("value", "variable"), c("bp_pred", "row_id"))

fake_predictions_m_spm <- validation_fake$predictions_m_spm
fake_predictions_m_long_spm <- melt(fake_predictions_m_spm, id.vars= c("id"))
setnames(fake_predictions_m_long_spm, c("value", "variable"), c("m_pred", "row_id"))

fake_predictions_long_spm <- as.data.table(merge(fake_predictions_bp_long_spm, fake_predictions_m_long_spm, by = c("row_id", "id")))
fake_predictions_long_spm[, model := "SPM"]



fake_predictions_bp_n <- validation_fake$predictions_bp_n
fake_predictions_bp_long_n <- melt(fake_predictions_bp_n, id.vars= "id")
setnames(fake_predictions_bp_long_n, c("value", "variable"), c("bp_pred", "row_id"))

fake_predictions_m_n<- validation_fake$predictions_m_n
fake_predictions_m_long_n <- melt(fake_predictions_m_n, id.vars= c("id"))
setnames(fake_predictions_m_long_n, c("value", "variable"), c("m_pred", "row_id"))

fake_predictions_long_n <- as.data.table(merge(fake_predictions_bp_long_n, fake_predictions_m_long_n, by = c("row_id", "id")))
fake_predictions_long_n[, model := "Naive"]

fake_predictions_long<-rbindlist(list(fake_predictions_long_spm, fake_predictions_long_n))


fake_predictions_median <- fake_predictions_long[, .(q500_bp = q500(bp_pred), 
                                                     mean_bp= mean(bp_pred),
                                    q500_m = q500(m_pred),
                                    mean_m = mean(m_pred)), keyby = .(id, model)]


convert_id_to_text <- function(id){
  return(glue::glue("Id = {id}"))
}
fake_predictions_median[, name_id := convert_id_to_text(id)]
fake_predictions_median



fake_predictions_long[, name_id:= convert_id_to_text(id)]

q <- ggplot(fake_predictions_long, aes(x = bp_pred, group = model, colour = model, linetype = model))
q <- q + geom_density()
q <- q + geom_vline(data = fake_predictions_median, aes(xintercept = mean_bp, colour = model, group = model, linetype = model))
q <- q + xlab(TeX("Precicted $BP_F$")) + ylab("Density") + theme(legend.title = element_blank(), legend.position = "bottom")
q + facet_wrap(vars(name_id), ncol = 1)


ggsave("images/validation_individual_examples_bp.pdf",
       width = 17,
       height = 10,
       units = "cm")

q <- ggplot(fake_predictions_long, aes(x = m_pred, group = model, colour = model, linetype = model))
q <- q + geom_density()
q <- q + geom_vline(data = fake_predictions_median, aes(xintercept = mean_m, colour = model, group = model))
q <- q + xlab(TeX("Precicted $p$")) + ylab("Density") + theme(legend.title = element_blank(), legend.position = "none")#, legend.position = "bottom")
#q <- q + scale_color_manual(values = "#F8766D")+xlim(0,1)
q + facet_wrap(vars(name_id), ncol = 1)




q <- ggplot(fake_predictions_long[model == "SPM"], aes(x = m_pred, group = model, colour = model, linetype = model))
q <- q + geom_density()
q <- q + geom_vline(data = fake_predictions_median[model == "SPM"], aes(xintercept = mean_m, colour = model, group = model))
q <- q + xlab(TeX("Precicted $p$")) + ylab("Density") + theme(legend.title = element_blank(), legend.position = "none")#, legend.position = "bottom")
q <- q + scale_color_manual(values = "#00BFC4") +xlim(0,1)
q + facet_wrap(vars(name_id), ncol = 1)
ggsave("images/validation_individual_examples_m_spm.pdf",
       width = 17,
       height = 10,
       units = "cm")

q <- ggplot(fake_predictions_long[model == "Naive"], aes(x = m_pred, group = model, colour = model, linetype = model))
q <- q + geom_density()
q <- q + geom_vline(data = fake_predictions_median[model == "Naive"], aes(xintercept = mean_m, colour = model, group = model))
q <- q + xlab(TeX("Precicted $p$")) + ylab("Density") + theme(legend.title = element_blank(), legend.position = "none")#, legend.position = "bottom")
q <- q + scale_color_manual(values = "#F8766D")+xlim(0,1)
q + facet_wrap(vars(name_id), ncol = 1)


ggsave("images/validation_individual_examples_m_naive.pdf",
       width = 17,
       height = 10,
       units = "cm")


###### Compare diff, reiduals, bias ########
nrow(data_val[missing==1])/nrow(data_val)
mean(predictions_agg_spm$diff_m)
mean(predictions_agg_n$diff_m)
mean(abs(na.omit(predictions_agg_spm$residuals_bp)))
mean(abs(na.omit(predictions_agg_n$residuals_bp)))

q <- ggplot(data= predictions_agg_spm, aes(x =q500_bp, y = residuals_bp, colour = "spm"))
q <- q + geom_point() 
q <- q + geom_point(data= predictions_agg_n, aes(x =q500_bp, y = residuals_bp, colour = "naive"))
q <- q + geom_smooth(method = "lm")
q <- q + geom_smooth(data= predictions_agg_n, aes(x =q500_bp, y = residuals_bp, colour = "naive"), method = "lm")
q

ggsave("images/validation_residuals_bp.pdf",
       width = 17,
       height = 10,
       units = "cm")

q <- ggplot(data= predictions_agg_spm, aes(x = diff_m, colour = "spm", fill = "spm"))
q <- q + geom_histogram(bins = 100, alpha = 0.5, aes(y = ..density..))
q <- q + geom_histogram(data= predictions_agg_n, aes(x = diff_m, y = ..density.., colour = "naive", fill = "naive"), bins = 100, alpha= 0.5)
q<- q + xlab("missing_status - mean probability of missing")
q

ggsave("images/validation_diff_missing_pred.pdf",
       width = 17,
       height = 10,
       units = "cm")

