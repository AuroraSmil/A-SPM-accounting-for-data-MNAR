library(INLA)
library(brinla)
#INLA:::inla.dynload.workaround()
library(tidyverse)
library(data.table)
library(latex2exp)
library(tikzDevice)

# Read data ----

data<- readRDS("/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")
# Set samples sample size: number of participants
sample_size <- 15000
index <- sample(seq(1,nrow(data)), sample_size, replace = FALSE )
data_small <- as.data.table(data[index])
n <- nrow(data_small)
data_sim <- data_small

# Set hyperpar for c 
c_mu = 0
c_sigma = 1
verbose = T
compute_all = FALSE


data_sim <- data_small
  
sigma2 = 0.001

#Prepare the response variables:
n <- nrow(data_sim)
y_gaussian <- c(data_sim$bp_3_corr,rep(NA,n))
m_binomial <- c(rep(NA,n),data_sim$missing)
joint_response <- list(y_gaussian,m_binomial)
  
#Specify number of trials for binomial repsonse variable
ntrials <- c(rep(NA,n),rep(1,n))

#Prepare covariates
covariates <- data.frame(alpha_0 = c(rep(1,n),rep(NA,n)),
                          beta_0 = c(rep(NA,n),rep(1,n)),
                          y_SEX = c(data_sim$sex,rep(NA,n)),
                          y_BP2 = c(data_sim$bp_2_corr,rep(NA,n)),
                          y_AGE = c(data_sim$age_2,rep(NA,n)),
                          y_BMI = c(data_sim$bmi_2,rep(NA,n)),
                          m_SEX = c(rep(NA,n),data_sim$sex),
                          m_BP2 = c(rep(NA,n),data_sim$bp_2_corr),
                          m_AGE = c(rep(NA,n),data_sim$age_2),
                          m_BMI = c(rep(NA,n),data_sim$bmi_2),
                          y_eps1 = c(1:n,rep(NA,n)),
                          m_eps1 = c(rep(NA,n),1:n)
)


joint_data <- c(covariates)

joint_data$Y <- joint_response

# SPM formula 
formula_spm = Y ~ -1 + alpha_0 + 
  y_SEX + 
  f(y_BP2,model="rw2",constr=T) +
  f(y_BMI,model="rw2",constr=T) +
  f(y_AGE,model="rw2",constr=T) +
  beta_0 + m_SEX + 
  f(m_BMI,model="rw2",constr=T) +
  f(m_AGE,model="rw2",constr=T) +
  f(m_BP2,model="rw2",constr=T) +
  f(y_eps1,model="iid") +
  f(m_eps1,copy="y_eps1",fixed=F,param=c(c_mu,c_sigma)) 
link= rep(NA, 2*(n ))
link[which(is.na(y_gaussian[1:(n )]))] = 1
link[(n) + which(is.na(m_binomial[((n)+1):(2*(n))]))] = 2
length(link)
start_time <- Sys.time()

#Run SPM 
if(compute_all == FALSE ){
  joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                     Ntrials=ntrials, verbose=verbose,
                     control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                     control.predictor = list(link = link),
                     control.compute = list(dic = TRUE))
}else{
  joint_inla <- inla(formula_spm, family=c("gaussian","binomial"),data=joint_data,
                     Ntrials=ntrials, verbose=verbose,
                     control.family=list(list(initial=log(1/sigma2^2),fixed=T),list()),
                     control.predictor = list(link = link),
                     control.compute=list(config=TRUE))
  
}

#Save model
saveRDS(joint_inla, file = "/home/aurorach/data_HUNT_aurora/master/non_lin_effects.RData")

# Some statistics
joint_inla$dic$dic
joint_inla$dic$p.eff
joint_inla$summary.fixed
summary(joint_inla)
bri.hyperpar.summary(joint_inla)

#Plot ---- 


joint_inla <- readRDS(file = "/home/aurorach/data_HUNT_aurora/master/non_lin_effects.RData")
data_random_age_missing <-  joint_inla$summary.random$m_AGE
data_random_age_bp <-  joint_inla$summary.random$y_AGE
data_random_bmi_missing <-  joint_inla$summary.random$m_BMI
data_random_bmi_bp <-  joint_inla$summary.random$y_BMI
data_random_BP_missing <-  joint_inla$summary.random$m_BP2
data_random_BP_bp <-  joint_inla$summary.random$y_BP2

override.linetype <- c(1,2,3,4,5,6)
q <- ggplot(data = data_random_age_missing, aes(x = ID))
q <- q + geom_line(aes(y = mean, colour = "age_missing", linetype = "age_missing"))
q <- q + geom_line(data = data_random_age_bp, aes(y = mean, colour = "age_bp", linetype = "age_bp"))
q <- q + geom_line(data = data_random_bmi_missing, aes(y = mean, colour = "bmi_missing", linetype = "bmi_missing"))
q <- q + geom_line(data = data_random_bmi_bp, aes(y = mean, colour = "bmi_bp", linetype = "bmi_bp"))
q <- q + geom_line(data = data_random_BP_missing, aes(y = mean, colour = "BP_missing", linetype = "BP_missing"))
q <- q + geom_line(data = data_random_BP_bp, aes(y = mean, colour = "BP_bp", linetype = "BP_bp"))
q <- q + scale_color_discrete(labels = unname(TeX(c("$\\f_{BP}(age)$",
                                                    "$\\f_{m}(age)$",
                                                    "$\\f_{BP}(BMI)$" ,
                                                    "$\\f_{m}(BMI)$" ,
                                                    "$\\f_{BP}(BP_c)$",
                                                    "$\\f_{m}(BP_c)$"))))
q <- q + labs(colour  = "Guide name", linetype = "Guide name")
q <- q + guides(colour = guide_legend(override.aes = list( linetype = override.linetype)))
q <- q + scale_linetype(guide = FALSE)
q <- q + xlab(element_blank())
q <- q + ylab("Effect") +guides(fill = FALSE) 
q <- q + xlim(-2.5, 2.5)
q <- q+ theme( legend.position = "bottom", legend.title = element_blank())


ggsave(glue::glue("images/non_lin_effects_participants_without_ci_{sample_size}.pdf"),
       q,
       width = 15,
       height = 20,
       units = "cm")

q <- ggplot(data = data_random_age_missing, aes(x = ID))
q <- q + geom_line(aes(y = mean, colour = "age_missing", linetype = "age_missing"))
q <- q + geom_ribbon(aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "age_missing"), alpha = 0.5)
q <- q + geom_line(data = data_random_age_bp, aes(y = mean, colour = "age_bp", linetype = "age_bp"))
q <- q + geom_ribbon(data = data_random_age_bp, aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "age_bp"), alpha = 0.5)
q <- q + geom_line(data = data_random_bmi_missing, aes(y = mean, colour = "bmi_missing", linetype = "bmi_missing"))
q <- q + geom_ribbon(data = data_random_bmi_missing, aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "bmi_missing"), alpha = 0.5)
q <- q + geom_line(data = data_random_bmi_bp, aes(y = mean, colour = "bmi_bp", linetype = "bmi_bp"))
q <- q + geom_ribbon(data = data_random_bmi_bp, aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "bmi_bp"), alpha = 0.5)
q <- q + geom_line(data = data_random_BP_missing, aes(y = mean, colour = "BP_missing", linetype = "BP_missing"))
q <- q + geom_ribbon(data = data_random_BP_missing, aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "BP_missing"), alpha = 0.5)
q <- q + geom_line(data = data_random_BP_bp, aes(y = mean, colour = "BP_bp", linetype = "BP_bp"))
q <- q + geom_ribbon(data = data_random_BP_bp, aes(ymin= `0.025quant`, ymax = `0.975quant`, fill = "BP_bp"), alpha = 0.5)
q <- q + scale_color_discrete(labels = unname(TeX(c("$\\f_{BP}(age)$",
                                                    "$\\f_{m}(age)$",
                                                    "$\\f_{BP}(BMI)$" ,
                                                    "$\\f_{m}(BMI)$" ,
                                                    "$\\f_{BP}(BP_c)$",
                                                    "$\\f_{m}(BP_c)$"))))
q <- q + labs(colour  = "Guide name", linetype = "Guide name")
q <- q + xlab(element_blank())+  xlim(-2.5, 2.5)
q <- q + ylab("Effect") +guides(fill = FALSE) 
q <- q + guides(colour = guide_legend(override.aes = list( linetype = override.linetype)))
q <- q + scale_linetype(guide = FALSE)
q <- q+ theme( legend.position = "bottom", legend.title = element_blank())

# Use contol.family to fix sigma2. Specified through log-precision.

ggsave(glue::glue("images/non_lin_effects_participants_with_ci_{sample_size}.pdf"),
       q,
       width = 15,
       height = 20,
       units = "cm")
end_time <- Sys.time()


