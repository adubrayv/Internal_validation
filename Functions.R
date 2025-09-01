####function
library(glmnet)
library(survival)
library(peperr)
library(pec)
library(caret)
library(rpart)
library(prodlim)
library(timeROC)

### Coefficients simulation
###100 positive and 100 negative coefficients using uniform distribution
simul_coefs = function(){
  set.seed(123)
  B1 = runif(n_positive_genes, 0.01,0.1)
  B2 = runif(n_negative_genes,-0.1,-0.01)
  B3 = rnorm(n_null_genes,0,0)
  coefficients = c(B1, B2, B3, 0.3, -0.5, 0.02, -0.8)
  return(coefficients)
}

###Data simulation
##mu, vv fixed in parameters_def
## first stage: fixation of log-TPM data sor each replicats
## second stage: fixation of clinical covariables: TNM/ HPV/ AGE/ SEX
## Calculation for each sample of each dataset using inverted Cox method a survival time
## Time point = 5 if survival time over 5 years censored = 0

simul_data = function(n,r) {
  set.seed(12 + r)
  i = which(n_sample == n)
  mu = abs(rnorm(n_genes,mean_mu,sd_mu))
  vv = rlnorm(n_genes, mean_vv ,sd_vv)
  dat = replicate(n,log(rsnorm(n_genes, mean = exp(mu), sd = exp(mu/scale_mu)*vv, xi)))
  dat[is.na(dat) | dat<0] = 0
  
  dat = as.data.frame(t(dat))
  dat$TNM = sample(1:4, n, replace = TRUE , prob = probs_TNM)
  dat$HPV = rbinom(n, 1, prob_HPV)
  dat$AGE = rnorm(n, mean_age, sd_age)
  dat$SEX = rbinom(n, 1, prob_sex)
  
  meanvector = sapply(dat, mean)
  H = -log(runif(n))  * exp(-(as.matrix(dat) %*% coefs) + sum(meanvector * coefs)) 
  dat$time = sapply(H, function(x) baseline_time[which(baseline_hazard > x)[1]])
  dat$time[is.na(dat$time)] = 5
  dat$event = 1*(dat$time<5)
}

##path
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##preparation for bootstrap and cross-validation,creation of survival object and optimization of lambda.1se by 10 fold cross-validation
survobj_fun = function(v){
  survival::Surv((v$surv), (v$event))
}
##calculation of lambda.1se
lambda_1se = function(v, survobj, a){
  fitlambda <- cv.glmnet(as.matrix(v[,1:15004]), survobj, family = "cox", type.measure = "C", nfolds = 10, alpha = a)
  bestlambda = fitlambda$lambda.1se
  return(bestlambda)
}

##best glmnet fit
fitmodel = function(v, a, survobj, bestlambda){
  fitmodel = glmnet(as.matrix(v[,1:15004]), alpha = a ,survobj, family = "cox", lambda = bestlambda)
  return(fitmodel)
}
                    

### oracle AUC
## selection of covariables associated with a simulated (naïve) non-zero coefficients
## Matrix of covariables values associated with non-zero coefficients
## calculation of linear predictor as test prediction for timeRoc package
oracle_fun <- function(v) {
  # Convert coefficients to matrix and subset non-zero coefficients
  matrix_coef <- as.data.frame(as.matrix(v[,1:15004]))
  NZindex <- which(coefs != 0)
  NZcoef <- coefs[NZindex]
  # Extract relevant columns for non-zero coefficients
  oracle_data <- v[, NZindex]
  # Compute the linear predictor (LP)
  LP <- as.matrix(oracle_data) %*% NZcoef
  v$LP = LP
  # Calculate time-dependent AUC using timeROC
  oracleAUC <- timeROC(
    T = v$time,
    delta = v$event,
    marker = LP,
    cause = 1,
    times = seq(1, 4.9, 0.1),
    ROC = TRUE
  )$AUC
  oracleAUC <- as.data.frame(oracleAUC)
  # Calculate the Integrated Brier Score using PEC
  pecfit <- pec(
    object = list(coxph(Surv(time, event) ~ LP, data = v, x = TRUE, y = TRUE)),
    formula = Surv(time, event) ~ LP,
    data = v,
    exact = FALSE,
    times = seq(0, 4.9, 0.1),
    cens.model = "cox",
    splitMethod = "none",
    B = 0
  )
  ibs <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0.1)
  # Return both AUC and IBS as a list
  return(list(oracleAUC = oracleAUC, IBS = unname(ibs[1, 1:50])))
}

##external validation
##model fitted using lambda.1se function for external validation
EV_AUC = function(model, bestlambda) {
  set.seed(1234)
  datatest = as.matrix(sim_valid[, 1:15004])
  sim_valid$preds = predict(model, newx = datatest, s = bestlambda)
  if (any(is.na(sim_valid$preds))) {
    next
  }
  AUCroc = timeROC(sim_valid$surv, sim_valid$event, sim_valid$preds, cause = 1, times = seq(1,4.99,0.1), ROC=TRUE)$AUC
  AUC_EV = as.data.frame(AUCroc)
  return(AUC_EV)
}

getAUC_EV <- function(AUC_EV_list) {
  AUC_EV_matrix <- do.call(cbind, AUC_EV_list)
  mean_AUC_EV <- rowMeans(AUC_EV_matrix)
  Q10_AUC_EV <- apply(AUC_EV_matrix, 1, quantile, probs = 0.1)
  Q90_AUC_EV <- apply(AUC_EV_matrix, 1, quantile, probs = 0.9)
  AUC_EV_df <- data.frame(mean_AUC_EV = mean_AUC_EV, Q10_AUC_EV = Q10_AUC_EV, Q90_AUC_EV = Q90_AUC_EV)
  return(AUC_EV_df)
}




