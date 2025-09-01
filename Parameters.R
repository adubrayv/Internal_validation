###Parameters of simulation
library(fGarch)
library(glmnet)
library(survival)
library(timeROC)

##Number of replication of scenario and sample size by scenario of simulation
reps = 100 ##number of replication of dataset per scenario
n_sample = c(25,50,75,100, 500, 1000) ##effective size of scenario 1/2/3
n = n_sample

##coefficients associated to the count data and the clinical data
##random effect uniform distribution for 100 transcripts
##100 negative coefficients

n_genes = 15000 ##number of transcript
n_positive_genes = 100 ## number of positive coefficients
n_negative_genes = 100 ## number of negative coefficients
n_null_genes = n_genes - n_positive_genes - n_negative_genes

##parameters of dispersion of count matrix for trancriptomic data using a log-normal asymetric distribution
mean_mu = 2.3 
sd_mu = 1.8
scale_mu = 1.3
mean_vv = log(2.8)
sd_vv = 0.4
xi = 2

##parameters of clinical data
probs_TNM = c(0.22,0.13,0.25,0.40) ## probability of stage I/II/III/IV of TNM staging
prob_HPV = 0.3 ##30% of HPV +
prob_sex = 0.3 ##30% of woman
mean_age = 65 ##mean age
sd_age = 10 ##standard deviation of age variable

##baseline hazard function of scandare study
load(file = "baseline_hazard_values.RData")
baseline_hazard <- baseline_hazard_values$hazard
baseline_time = baseline_hazard_values$time
