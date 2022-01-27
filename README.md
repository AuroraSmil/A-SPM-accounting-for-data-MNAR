# A shared parameter model accounting for dropout not at random in a predictive model for systolic blood pressure using data from The HUNT Study

## Author Contributions Checklist Form
<!--HOW TO COMPLETE THIS FORM:-->

<!--
1. Checkboxes in this document appear as follows: 

- [ ] This is a checkbox 

To check a checkbox, replace [ ] by [x], as follows: 

- [x] This is a checked checkbox 

Note that older versions of RStudio (versions lower than 1.3) may not create a formatted checkbox but will leave the original characters, i.e., literally "[ ]" or "[x]". It's fine to submit a PDF in this form.
 
2. For text answers, simply type the relevant text in the areas indicated. A blank line starts a new paragraph. 
 
3. Comments (like these instructions) provide additional instructions throughout the form. There is no need to remove them; they will not appear in the compiled document. 

4. If you are comfortable with Markdown syntax, you may choose to include any Markdown-compliant formatting in the form. For example, you may wish to include R code chunks and compile this document in R Markdown.
-->


### Data

#### Abstract

<!--
Provide a short (< 100 words), high-level description of the data
-->

This work focuses on a case study with a predictive model for systolic blood pressure (BP) using data from a longitudinal population-based health survey, the Trøndelag Health (HUNT) Study. Elevated BP increases the risk of developing diseases related to the brain, heart, blood vessels, and kidney. This medical condition affects more than 1.1 billion people. It accounts for over 10.8 million deaths per year, surpassing smoking as the leading preventable cause of death for middle-aged and older adults worldwide. Early detection, prevention, and treatment of elevated BP are of high priority in public health strategies. Thus, obtaining unbiased, accurate models for predicting future BP is of great interest in medical research. Indeed, this work is motivated by longitudinal population-based health surveys in medical research and clinical observational studies with a follow-up design and is inspired by a multidisciplinary collaboration on the epidemiology and etiology of elevated BP using data from the HUNT Study.


#### Availability

Data **cannot be made** publicly available.


<!-- If data are available by request to the authors or some other data owner, please make sure to explain the process of requesting access to the data. -->


<!--
The Journal of the American Statistical Association requires authors to make data accompanying their papers available to the scientific community except in cases where: 1) public sharing of data would be impossible, 2) suitable synthetic data are provided which allow the main analyses to be replicated (recognizing that results may differ from the "real" data analyses), and 3) the scientific value of the results and methods outweigh the lack of reproducibility.

Please discuss the lack of publicly available data. For example:
-	why data sharing is not possible,
-	what synthetic data are provided, and 
-	why the value of the paper's scientific contribution outweighs the lack of reproducibility.
-->

The HUNT data contains sensitive information and is protected by law. The following paragraph is from HUNTs guidelines for publication https://www.ntnu.edu/documents/140075/1295406997/Guidelines+for+publication+of+research+results+using+HUNT-data.pdf/007a5d98-4369-94a6-06b7-6264367b3faa?t=1600862719853:

"The Trøndelag Health Study (HUNT) has invited persons aged 13 - 100 years to four surveys between 1984
and 2019. Comprehensive data from more than 140,000 persons having participated at least once and
biological material from 78,000 persons are collected. The data are stored in HUNT databank and biological
material in HUNT biobank. HUNT Research Centre has permission from the Norwegian Data Inspectorate
to store and handle these data. The key identification in the database is the personal identification number
given to all Norwegians at birth or immigration, whilst de-identified data are sent to researchers upon
approval of a research protocol by the Regional Ethical Committee and HUNT Research Centre. To protect
participants' privacy, HUNT Research Centre aims to limit storage of data outside HUNT databank, and
cannot deposit data in open repositories. HUNT databank has precise information on all data exported to
different projects and are able to reproduce these on request. There are no restrictions regarding data export
given approval of applications to HUNT Research Centre. For more information see:
http://www.ntnu.edu/hunt/data"

Therefore we provide synthetic data with similar properties containing the same variables. This data ensures that all code can be run and inspected. 
However, the results will not be identical to those produced on the actual data. 

Still, we find it of scientific interest to explore how the models used in this work on real data can improve our understanding of blood pressure and make good predictive models for this condition. 

#### Description


The HUNT Study is a longitudinal population-based health survey in central Norway. Every adult citizen in the now-former county of Nord-Trøndelag was invited to participate in the first survey in 1984-86 (HUNT1). Since then, all adult residents in the screening area have been invited to clinical examinations and questionnaires in 1995-97 (HUNT2), 2006-08 (HUNT3), and 2017-19 (HUNT4).
This work uses HUNT data in a predictive model of BP eleven years ahead (BP_F) based on current BP (BP_I), age, sex, and BMI. 
As large proportions (33% -43%) of participants are lost to follow-up in between consecutive surveys, we want to account for data being missing not at random (MNAR). Hence we construct the variable m indicating if a participant is missing (m=1) or present (m= 0). 
The dataset contains the following variables: BP_F, m, sex, age, BMI, and BP_I.

In this work, we define two cohorts for which we only consider participants with complete records of the explanatory variables. The modeling cohort (HUNT2 cohort, n = 64 386) consists of explanatory variables from HUNT2 and response values from HUNT3. 43.1% of the participants drop out prior to HUNT3. The validation cohort (HUNT3 cohort, n = 50 807) consists of explanatory variables from HUNT3 and response variables from HUNT4. In the HUNT3 cohort, 33.3% drop out prior to HUNT4. When needed, we use a subscript to indicate which HUNT survey the variable originates from (e.g., BP_2 denotes BP at HUNT2).

The simulated dataset for reproducibility contains the same set of variables and is constructed based on the HUNT2 cohort. 


#### Additional Information (optional)
<!-- 
OPTIONAL: Provide any additional details that would be helpful in understanding the data. If relevant, please provide unique identifier/DOI/version information and/or license/terms of use.
-->
We have created two simulated datasets, the first of the same size as the original data and the second consisting only of 1000 simulated participants. To run the code on a dataset of the original size, one needs a lot of computational power. We had 32 CPUs and 32 G RAM available. However, the small dataset can be run on a desktop laptop to inspect the code. We note that the small dataset will not provide even remotely similar results as the actual data because the models used need datasets of certain sizes to establish the proper connections. 


###  Code

#### Abstract

<!--
Provide a short (< 100 words), high-level description of the code. If necessary, more details can be provided in files that accompany the code. If no code is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

We use integrated nested Laplace approximations (INLA) to fit both a shared parameter model (SPM) and a naive model. INLA supports fitting models with multiple likelihoods, which is needed when fitting the SPM. This method falls within the Bayesian framework. In contrast to regular MCMC, INLA uses numerical approximations instead of samples. 

For both the SPM and the naive model, we use normal priors for the latent field and gamma priors for the hyperparameters. 

We sample from the posterior distributions of all the parameter estimates for model predictions and use these to obtain the posterior predictive distributions for all participants. 

All code needed to reproduce the work done is available at https://github.com/AuroraSmil/A-SPM-accounting-for-data-MNAR. 



#### Description

Fitting and sampling from the models used in this work are computationally demanding. Therefore we have split the different parts of this work into different scripts. 

##### Fitting Models Using INLA
To fit the SPM and naive models, run code/INLA/inla_spm_naive.R. Make sure to edit the storage location at the beginning of the script.
The results are compared, and the effect of the association parameter is explored at the end of the script. The comparison and exploration part can be run separately if the packages at the beginning of the script are loaded.  

We also added code to explore modeling all continuous explanatory models as additive effects in the SPM. To do so, run code/INLA/non_lin_effects_spm.R. This is extremely computationally demanding. Hence, we recommend starting with the smaller dataset. 

##### Validation of Model Predictions and MNAR Assumption

code/validation/validation_func.R contains functions to obtain model predictions, computing CRPS and Brier scores, and obtaining model predictions given missing status. 

To obtain model predictions and all plots and tables related to this, run code/validation/validation_predictions.R.
To obtain the results related to the validation of the MNAR assumption, run code/validation/validation_MNAR.R.

##### Simulatoin Studies

We perform several simulation studies for which we construct new values for BP_F and m based on given parameter values. General functions can be found in code/simulation_study/model_functions.R. This script contains general functions for fitting the SPM, naive BP model, and naive missing model. In addition, it contains a function to generate simulated data based on given parameter estimates and a function for plotting the results of the simulation studies.

To explore how the models perform on data simulated with known parameters, we performed several simulation studies with posterior mean estimates from the SPM and the naive model with data MAR and MNAR. The code for this can be found in code/simulation_study/simulation_study_reprodicability.R
The figures are created at the end of the script and can also be obtained by running code/simulation_study/read_simulation_results_reproduacbility.R. 

We also perform a simulation study to explore the two validation schemes proposed in this work.
We explore how the CRPS and Brier score varies between present and missing data. This is done by simulating a dataset with explanatory variables from the HUNT2 cohort before refitting both the SPM and naive model to this dataset. This is done in code/simulation_study/refit_models_for_validation_spm_parameters.R and code/simulation_study/refit_models_for_validation_naive_parameters.R. 
Then these models are used to predict BP_F and m on simulated data with explanatory variables from the HUNT3 cohort and compute the CRPS and Brier score for all participants. This is done in code/simulation_study/simulation_study_validation_predictions. To read the results, run code/simulation_study/read_simulatoin_results_validation_predctions.R followed by code/simulation_study/evaluate_validation_prediction.R. 

We also perform a simulation study to explore the differences we can expect when predicting BP_f given missing status and without this knowledge. To reproduce this run code/simulation_study/simulation_study_mnar_spm.R and code/simulation_study/simulation_study_mnar_naive.R to explore both the case when the data is MNAR and MAR respectively. Evaluate the results by running code/simulation_study/evaluate_validatoin_mnar.R. code/simulation_study/create_simulated_HUNT3_cohort_data_for_validation.R creates n simulated datasets. 

##### Prior Sensitivity Analysis

We performed a small prior sensitivity analysis to explore the robustness concerning prior misspecification of the association parameter c. To reproduce this run code/prior_sisitivity/prior_sensitivity_c.R.

##### Preprocessing 
The folder data_clean contains the code for preprocessing the data. 

#### Version of primary software used

R version 4.0.4 (2021-02-15)

#### Libraries and dependencies used by the code


<!--
Include version numbers (e.g., version numbers for any R or Python packages used)
-->

attached base packages:
stats4    parallel  grid      stats     graphics  grDevices utils     datasets 
methods   base

other attached packages:
corrplot_0.92       fGarch_3042.83.2    fBasics_3042.89.1   timeSeries_3062.100
timeDate_3043.102   sn_1.6-2            Hmisc_4.5-0         Formula_1.2-4      
survival_3.2-7      lattice_0.20-41     lubridate_1.7.9.2   mgcv_1.8-33        
nlme_3.1-152        scoringRules_1.0.1  DescTools_0.99.41   glue_1.4.2         
forcats_0.5.1       stringr_1.4.0       dplyr_1.0.4         purrr_0.3.4        
readr_1.4.0         tidyr_1.1.2         tibble_3.0.6        tidyverse_1.3.0    
brinla_0.1.0        INLA_21.02.23       sp_1.4-5            foreach_1.5.1      
Matrix_1.3-2        gridtext_0.1.4      tikzDevice_0.12.3.1 latex2exp_0.4.0    
ggplot2_3.3.3       gridExtra_2.3       cowplot_1.1.1       foreign_0.8-81     
data.table_1.13.6  


#### Supporting system/hardware requirements (optional)

<!--
OPTIONAL: System/hardware requirements including operating system with version number, access to cluster, GPUs, etc.
-->

Ubuntu 18.04.6 LTS


### Additional information (optional)

<!--
OPTIONAL: By default, submitted code will be published on the JASA GitHub repository (http://github.com/JASA-ACS) as well as in the supplementary material. Authors are encouraged to also make their code available in a public code repository, such as on GitHub, GitLab, or BitBucket. If relevant, please provide unique identifier/DOI/version information (e.g., a Git commit ID, branch, release, or tag). If the code and workflow are provided together, this section may be omitted, with information provided in the "Location" section below.
-->

###  Reproducibility workflow

<!--
The materials provided should provide a straightforward way for reviewers and readers to reproduce analyses with as few steps as possible. 
-->

### Workflow

#### Instructions

Load all necessary packages. Check especially that INLA is loaded correctly. 

To reproduce the graphs and numbers for all tables in this work, follow the above order and run all scripts. 
To run the code on a similar size as this work, one needs a lot of memory and the possibility for parallelization. Even then, it is a time-consuming process. 

If the aim is to check the code and thought processes behind the modeling, we recommend using a much smaller dataset. Either way, expect runtimes over 8 hours to reproduce the entire work. 

<!--
Describe how to use the materials provided to reproduce analyses in the manuscript. Additional details can be provided in file(s) accompanying the reproducibility materials. If no workflow is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

<!--
### Additional information (optional)
-->
<!--
OPTIONAL: Additional documentation provided (e.g., R package vignettes, demos or other examples) that show how to use the provided code/software in other settings.
-->

# Notes (optional)

<!--
OPTIONAL: Any other relevant information not covered on this form. If reproducibility materials are not publicly available at the time of submission, please provide information here on how the reviewers can view the materials.
-->
