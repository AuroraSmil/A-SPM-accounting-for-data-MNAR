library(data.table)
library(foreign)

# Read raw data ----
data_exl_age<- read.spss("/home/aurorach/data_HUNT_aurora/master/2021-11-16_111751_Data.sav", 
                         to.data.frame = TRUE,
                         use.value.labels = FALSE)

data_exl_age <- as.data.table(data_exl_age)



data_age<- read.spss("/home/aurorach/data_HUNT_aurora/master/2021-11-17_111751_Tillegg.sav", 
                     to.data.frame = TRUE,
                     use.value.labels = FALSE)

data_age <- as.data.table(data_age)

# Change names ----
data_exl_age_interest<- data_exl_age[, .(`PID.111751`,
                          Sex,
                          `BPSystMn12.NT1BLM`,
                          `BPMedEv.NT1BLQ1`,
                          `BPMedCu.NT1CvdQ`,
                          `Bmi.NT1BLM`,
                          `BPSystMn23.NT2BLM`,
                          `BPMedCu.NT2BLQ1`,
                          `Bmi.NT2BLM`,
                         `RegisStat`,
                         `RegisStatDat`,
                         `ObsEndDat`,
                         `BirthYear`,
                         `BPSystMn23.NT3BLM`, 
                         `BPMedEv.NT3BLQ1`,
                         `BPMedSiEffEv.NT3CvdQ`,
                         `Bmi.NT3BLM`,
                         `BPSystMn23.NT4BLM`,
                         `BPMedCu.NT4BLQ1`,
                         `MedPrescCu.NT4BLQ1`)]

setnames(data_exl_age_interest,
         c("PID.111751",
           "RegisStat",
           "RegisStatDat",
           "ObsEndDat",
           "BirthYear",
           "Sex",
           "Bmi.NT1BLM",
           "BPSystMn12.NT1BLM",
           "BPMedEv.NT1BLQ1",
           "BPMedCu.NT1CvdQ",
           "Bmi.NT2BLM",
           "BPSystMn23.NT2BLM", 
           "BPMedCu.NT2BLQ1",
           "Bmi.NT3BLM"  ,
           "BPSystMn23.NT3BLM",
           "BPMedSiEffEv.NT3CvdQ",
           "BPMedEv.NT3BLQ1",
           "BPSystMn23.NT4BLM",
           "MedPrescCu.NT4BLQ1",
           "BPMedCu.NT4BLQ1"
           ),
         c("id",
           "register_status",
           "register_status_date",
           "observed_end_date",
           "birth_year",
           "sex",
           "bmi_1",
           "bp_1", 
           "bp_use_of_med_1", 
           "bp_use_med_now_1",
           "bmi_2",
           "bp_2",
           "bp_use_med_now_2",
           "bmi_3",
           "bp_3",
           "bp_pain_from_use_med_now_3",
           "bp_use_of_med_3",
           "bp_4",
           "use_of_any_med_4",
           "bp_use_med_now_4"
           ))



data_age_interest <- data_age[, .(`PID.111751`,
                              `PartAg.NT1BLQ1`,
                              `PartAg.NT2BLQ1`,
                              `PartAg.NT3BLQ1`,
                              `PartAg.NT4BLQ1`)]


setnames(data_age_interest,
         c("PID.111751",
           "PartAg.NT1BLQ1",
           "PartAg.NT2BLQ1",
           "PartAg.NT3BLQ1",
           "PartAg.NT4BLQ1"
         ),
         c("id",
           "age_1",
           "age_2",
           "age_3",
           "age_4"
         ))

### Merge data -----
data_interest <- merge(data_exl_age_interest, data_age_interest, on = "id")

# Check no data is lost

nrow(data_exl_age_interest)
nrow(data_age_interest)
nrow(data_interest)

# Converts dates to our calender 
date_function <- function(x){
  return(as.Date(x/86400, origin = "1582-10-14"))
} 

# Check dead and moved ----
date_start_hunt3 <- as.Date("2006-10-03")
date_start_hunt4 <- as.Date("2017-09-01")
data_interest[, register_status_date_convert := date_function(register_status_date)]
data_interest[, observed_end_date_convert := date_function(observed_end_date)]
data_interest[, dead_3 := 0]
data_interest[register_status== 5 & register_status_date_convert< date_start_hunt3, dead_3 := 1]
data_interest[, moved_3 := 0]
data_interest[register_status %in% c(1,2,3) & !is.na(observed_end_date) & observed_end_date_convert < date_start_hunt3, moved_3 := 1]
data_interest[, dead_4 := 0]
data_interest[register_status== 5 & register_status_date_convert< date_start_hunt4, dead_4 := 1]
data_interest[, moved_4 := 0]
data_interest[register_status %in% c(1,2,3) & !is.na(observed_end_date) & observed_end_date_convert < date_start_hunt4, moved_4 := 1]

#People without full records from HUNT2
nrow(data_interest[!is.na(age_2) & !is.na(sex) & ( is.na(bp_2) | is.na(bmi_2))])
nrow(data_interest[!is.na(age_2) & !is.na(sex) & !is.na(bp_2) & is.na(bmi_2)])
nrow(data_interest[!is.na(age_3) & !is.na(sex) & ( is.na(bp_3) | is.na(bmi_3))])


nrow(data_interest[is.na(register_status),])
nrow(data_interest[dead_3 == 1])
nrow(data_interest[moved_3 == 1])
nrow(data_interest[dead_4 == 1])
nrow(data_interest[moved_4 == 1])

#Correct for BP medication ----

data_interest[, bp_2_corr:= bp_2]
data_interest[, bp_3_corr:= bp_3]
data_interest[, bp_4_corr:= bp_4]
data_interest[bp_use_med_now_2 == 1, bp_2_corr := bp_2 + 15]
data_interest[!is.na(bp_pain_from_use_med_now_3), bp_3_corr := bp_3 +15] #Not really correct but best we got.. 
data_interest[use_of_any_med_4 == 1 &bp_use_med_now_4 ==1, bp_4_corr := bp_4+ 15]

nrow(data_interest[bp_use_med_now_2 == 1])
nrow(data_interest[!is.na(bp_pain_from_use_med_now_3)])
nrow(data_interest[use_of_any_med_4 == 1 &bp_use_med_now_4 ==1])

# Create HUNT23 data----
data_H23 <- copy(na.omit(data_interest, cols = c("age_2", "sex", "bp_2_corr", "bmi_2")))
nrow(data_H23)


nrow(data_interest[bp_use_med_now_2 == 1])/nrow(data_H23)
#Create drop-out variable ----
#(1 if missing, 0 if present)

data_H23[, missing:= 0]
data_H23[is.na(bp_3_corr), missing := 1]

psych::describeBy(data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)])
psych::describeBy(data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)], 
                  data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)]$missing)
psych::describeBy(data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)][sex == 0], 
                  data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)][sex == 0]$missing)
psych::describeBy(data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)][sex == 1], 
                  data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)][sex == 1]$missing)

nrow(data_H23[missing == 1,])
nrow(data_H23[dead_3 == 1,])
nrow(data_H23[moved_3 == 1,])
nrow(data_H23[(dead_3 == 1 | moved_3 == 1) & missing == 0,])
nrow(data_H23[(dead_3 == 0 & moved_3 == 0) & missing == 1,])

nrow(data_H23[moved_3 == 1,]) + 
  nrow(data_H23[dead_3 == 1, ]) + 
  nrow(data_H23[(dead_3 == 0 & moved_3 == 0) & missing == 1,]) + 
  nrow(data_H23[is.na(register_status),])

saveRDS(data_H23,file="/home/aurorach/data_HUNT_aurora/master/HUNT23_adjusted_meds.RData")

# Standardize the data
data_H23_scaled <- copy(data_H23[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)])
data_H23_scaled[,bp_2_corr:= scale(bp_2_corr)]
data_H23_scaled[,age_2:= scale(age_2)]
data_H23_scaled[,bmi_2:= scale(bmi_2)]
data_H23_scaled[,bp_3_corr:= scale(bp_3_corr)]

psych::describeBy(data_H23_scaled)
saveRDS(data_H23_scaled,file="/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled.RData")


# True missing #####


data_H23_true_missing <- copy(data_H23[dead_3 != 1 & moved_3 != 1,])
nrow(data_H23_true_missing)
nrow(data_H23) -nrow(data_H23_true_missing)
### SUMMARY STATISTICS #####

grouped_summary_val <- psych::describeBy(data_H23_true_missing[, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing)], data_H23_true_missing$missing)
summary_val <- psych::describeBy(data_H23_true_missing[, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing)])

psych::describeBy(data_H23_true_missing[, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing , dead_3, moved_3)], data_H23_true_missing$sex)
psych::describeBy(data_H23_true_missing[sex == 0, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing , dead_3, moved_3)], data_H23_true_missing[sex == 0]$missing)
psych::describeBy(data_H23_true_missing[sex == 1, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing , dead_3, moved_3)], data_H23_true_missing[sex == 1]$missing)

# Scaled the data
data_H23_true_missing_scaled <- copy(data_H23_true_missing[, .(id, bp_2_corr, bp_3_corr, age_2, bmi_2, sex, missing)])
data_H23_true_missing_scaled[,bp_2_corr:= scale(bp_2_corr)]
data_H23_true_missing_scaled[,age_2:= scale(age_2)]
data_H23_true_missing_scaled[,bmi_2:= scale(bmi_2)]
data_H23_true_missing_scaled[,bp_3_corr:= scale(bp_3_corr)]

psych::describeBy(data_H23_true_missing_scaled)
saveRDS(data_H23_true_missing_scaled,file="/home/aurorach/data_HUNT_aurora/master/HUNT23_scaled_true_missing.RData")


#Create hunt34 data -----
data_H34 <- copy(na.omit(data_interest, cols = c("age_3", "sex", "bp_3_corr", "bmi_3")))
nrow(data_H34)


#Create drop-out variable ----
#(1 if missing, 0 if not)

data_H34[, missing:= 0]
data_H34[is.na(bp_4_corr), missing := 1]

psych::describeBy(data_H34[, .( bp_4_corr, bp_3_corr,  age_3,  bmi_3,sex, missing)])
psych::describeBy(data_H34[, .(bp_4_corr, bp_3_corr,  age_3,  bmi_3,sex, missing)], data_H34$sex)
psych::describeBy(data_H34[, .( bp_4_corr, bp_3_corr,  age_3,  bmi_3,sex, missing)], data_H34$missing)
psych::describeBy(data_H34[sex == 0, .( bp_4_corr, bp_3_corr,  age_3,  bmi_3,sex, missing)], data_H34[sex == 0]$missing)
psych::describeBy(data_H34[sex == 1, .( bp_4_corr, bp_3_corr,  age_3,  bmi_3,sex, missing)], data_H34[sex == 1]$missing)


nrow(data_H34[missing == 1,])
nrow(data_H34[dead_4 == 1,])
nrow(data_H34[moved_4 == 1,])
nrow(data_H34[(dead_4 == 1 | moved_4 == 1) & missing == 0,])
nrow(data_H34[(dead_4 == 0 & moved_4 == 0) & missing == 1,])

nrow(data_H34[moved_4 == 1,]) + 
  nrow(data_H34[dead_4 == 1, ]) + 
  nrow(data_H34[(dead_4 == 0 & moved_4 == 0) & missing == 1,]) + 
  nrow(data_H34[is.na(register_status),])

saveRDS(data_H34,file="/home/aurorach/data_HUNT_aurora/master/HUNT34_adjusted_meds.RData")

# Standardize the data with HUNT23 scale ----
data_H34_scaled <- copy(data_H34[, .(id, bp_4_corr, bp_3_corr, age_3, bmi_3, sex, missing)])
data_H34_scaled[,bp_3_corr:= scale(bp_3_corr,
                                  center = attributes(data_H23_scaled$bp_2_corr)$`scaled:center`[[1]], 
                                  scale = attributes(data_H23_scaled$bp_2_corr)$`scaled:scale`[[1]])]
data_H34_scaled[,age_3:= scale(age_3,
                               center = attributes(data_H23_scaled$age_2)$`scaled:center`[[1]], 
                               scale = attributes(data_H23_scaled$age_2)$`scaled:scale`[[1]])]
data_H34_scaled[,bmi_3:= scale(bmi_3,
                               center = attributes(data_H23_scaled$bmi_2)$`scaled:center`[[1]], 
                               scale = attributes(data_H23_scaled$bmi_2)$`scaled:scale`[[1]])]
data_H34_scaled[,bp_4_corr:= scale(bp_4_corr,
                                   center = attributes(data_H23_scaled$bp_3_corr)$`scaled:center`[[1]], 
                                   scale = attributes(data_H23_scaled$bp_3_corr)$`scaled:scale`[[1]])]

psych::describeBy(data_H34_scaled)
saveRDS(data_H34_scaled,file="/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23.RData")


# HUNT34 true mising----
data_H34_true_missing <- copy(data_H34[dead_4 != 1 & moved_4 != 1,])
nrow(data_H34_true_missing)
nrow(data_H34) -nrow(data_H34_true_missing)
### SUMMARY STATISTICS #####

grouped_summary_val <- psych::describeBy(data_H34_true_missing[, .( sex, age_3, bp_3_corr, bmi_3, bp_4_corr, missing)], data_H34_true_missing$missing)
summary_val <- psych::describeBy(data_H34_true_missing[, .( sex, age_3, bp_3_corr, bmi_3, bp_4_corr, missing)])

psych::describeBy(data_H34_true_missing[, .( sex, age_3, bp_3_corr, bmi_3, bp_4_corr, missing , dead_4, moved_4)], data_H34_true_missing$sex)
psych::describeBy(data_H34_true_missing[sex == 0, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing , dead_4, moved_4)], data_H34_true_missing[sex == 0]$missing)
psych::describeBy(data_H34_true_missing[sex == 1, .( sex, age_2, bp_2_corr, bmi_2, bp_3_corr, missing , dead_4, moved_4)], data_H34_true_missing[sex == 1]$missing)

# Standardize the data with HUNT23 values ----
data_H34_true_missing_scaled <- copy(data_H34_true_missing[, .(id, bp_3_corr, bp_4_corr, age_3, bmi_3, sex, missing)])
data_H34_true_missing_scaled[,bp_3_corr:= scale(bp_3_corr,
                                                center = attributes(data_H23_true_missing_scaled$bp_2_corr)$`scaled:center`[[1]], 
                                                scale = attributes(data_H23_true_missing_scaled$bp_2_corr)$`scaled:scale`[[1]])]
data_H34_true_missing_scaled[,age_3:= scale(age_3,
                                            center = attributes(data_H23_true_missing_scaled$age_2)$`scaled:center`[[1]], 
                                            scale = attributes(data_H23_true_missing_scaled$age_2)$`scaled:scale`[[1]])]
data_H34_true_missing_scaled[,bmi_3:= scale(bmi_3,
                                            center = attributes(data_H23_true_missing_scaled$bmi_2)$`scaled:center`[[1]], 
                                            scale = attributes(data_H23_true_missing_scaled$bmi_2)$`scaled:scale`[[1]])]
data_H34_true_missing_scaled[,bp_4_corr:= scale(bp_4_corr,
                                                center = attributes(data_H23_true_missing_scaled$bp_3_corr)$`scaled:center`[[1]], 
                                                scale = attributes(data_H23_true_missing_scaled$bp_3_corr)$`scaled:scale`[[1]])]


psych::describeBy(data_H34_true_missing_scaled)
saveRDS(data_H34_true_missing_scaled,file="/home/aurorach/data_HUNT_aurora/master/HUNT34_scaled_H23_true_missing.RData")


