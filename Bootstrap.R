## source
source("lambda1setuning.R")
##library
library(parallel)
library(foreach)
library(doParallel)
## Settings of bootstrap: original B = 100, a = 0.95/0.5/0.05
### 1) Original bootstrap using lambda.1se tuned in model selection for each scenario and penalized method
### 2) Extraction of the mean of AUC(t), Q10, Q90

##scenario1###
  ##lasso-like
  BT_lasso1 = lapply(1:100, function(i) AUC_bootstrap(sim_data100[[i]], 
                                                      0.95, 
                                                      bestlambdaenet1[[i]], 
                                                      15000, 
                                                      100, 
                                                      Weight = FALSE))
  
  AUC_BT_lasso1 = Get_AUC_BT(BT_lasso1)
  save(AUC_BT_lasso1, file = "result_BT_lasso1.RData")
  
  ##enet
  BT_enet1 <- lapply(1:100, function(i) AUC_bootstrap(sim_data100[[i]], 
                                                      0.5, 
                                                      bestlambdaenet1[[i]], 
                                                      15000, 
                                                      100, 
                                                      Weight = FALSE))
  

  AUC_BT_enet1 = Get_AUC_BT(BT_enet1)
  save(AUC_BT_enet1, file = "result_BT_enet1.RData")
  
  ##ridge-like
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data100", "bestlambdaridge1", "AUC_bootstrap"))
  BTpar_ridge1 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data100[[i]], 
                  0.05, 
                  bestlambdaridge1[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_ridge1 = Get_AUC_BT(BTpar_ridge1)
  save(AUC_BT_ridge1, file = "result_BT_ridge1.RData")
  
  ##adaptative Enet
  AUC_BT_adapt1 = lapply(1:100, function (i) AUC_bootstrap(sim_data100[[i]], 
                                                           0.5, 
                                                           bestlambdaenet1[[i]], 
                                                           15000, 
                                                           100, 
                                                           Weight = TRUE))
  result_BT_adapt1 = Get_AUC_BT(AUC_BT_adapt1)
  save(result_BT_adapt1, file = "result_BT_adapt1.RData")
  
  result_BT1 = cbind(AUC_BT_lasso1, AUC_BT_enet1, AUC_BT_ridge1)
  save(result_BT1, file="result_BT1.RData")
  
###scenario2##
    ##lasso-like 
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data500", "bestlambda2", "AUC_bootstrap"))
  BT_lasso2 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data500[[i]], 
                  0.95, 
                  bestlambda2[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_lasso2 = Get_AUC_BT(BT_lasso2)
  save(AUC_BT_lasso2, file = "result_BT_lasso2.RData")
  
  ##enet
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data500", "bestlambdaenet2", "AUC_bootstrap"))
  BT_enet2 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data500[[i]], 
                  0.5, 
                  bestlambdaenet2[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_enet2 = Get_AUC_BT(BT_enet2)
  save(AUC_BT_enet2, file = "AUC_BT_enet2.RData")
  
  ##ridge-like
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data500", "bestlambdaridge2", "AUC_bootstrap"))
  BT_ridge2 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data500[[i]], 
                  0.05, 
                  bestlambdaridge2[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_ridge2 = Get_AUC_BT(BT_ridge2)
  save(AUC_BT_ridge2, file = "AUC_BT_ridge2.RData")
  
  ##adatpative Enet
  AUC_BT_adapt2 = lapply(1:100, function (i) split_sample_AUC(sim_data500[[i]], 
                                                              0.5, 
                                                              bestlambdaenet2[[i]], 
                                                              15000, 
                                                              Weight = TRUE))
  result_BT_adapt2 = Get_AUC_BT(AUC_BT_adapt2)
  save(result_BT_adapt2, file = "result_BT_adapt2.RData")
  
  result_BT2 = cbind(AUC_BT_lasso2, AUC_BT_enet2, AUC_BT_ridge2)
  save(result_BT2, file="result_BT2.RData")

###scenario3##
    ##lasso-like
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data1000", "bestlambda3", "AUC_bootstrap"))
  BT_lasso3 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data1000[[i]], 
                  0.95, 
                  bestlambda3[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_lasso3 = Get_AUC_BT(BT_lasso3)
  save(AUC_BT_lasso3, file = "result_BT_lasso3.RData")
  
  ##enet
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data1000", "bestlambdaenet3", "AUC_bootstrap"))
  BT_enet3 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data1000[[i]], 
                  0.5, 
                  bestlambdaenet3[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_enet3 = Get_AUC_BT(BT_enet3)
  save(AUC_BT_enet3, file = "AUC_BT_enet3.RData")
  
  ##ridge-like
  numCores <- 4
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(survival)
    library(timeROC)
    library(riskRegression)
    library(pec)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data1000", "bestlambdaridge3", "AUC_bootstrap"))
  BT_ridge3 <- parLapply(cl, 1:100, function(i) {
    AUC_bootstrap(sim_data1000[[i]], 
                  0.05, 
                  bestlambdaridge3[[i]], 
                  15000, 
                  100, 
                  Weight = FALSE)
  })
  stopCluster(cl)
  AUC_BT_ridge3 = Get_AUC_BT(BT_ridge3)
  save(AUC_BT_ridge3, file = "AUC_BT_ridge3.RData")
  
  ##adaptative Enet
  AUC_BT_adapt3 = lapply(1:100, function (i) split_sample_AUC(sim_data1000[[i]], 
                                                              0.5, 
                                                              bestlambdaenet3[[i]], 
                                                              15000, 
                                                              Weight = TRUE))
  result_BT_adapt3 = Get_AUC_BT(AUC_BT_adapt3)
  save(result_BT_adapt3, file = "result_BT_adapt3.RData")
  
  result_BT3 = cbind(AUC_BT_lasso3, AUC_BT_enet3, AUC_BT_ridge3)
  save(result_BT3, file="result_BT3.RData")

    ## save the results
  result_BT1 = cbind(AUC_BT_lasso1, AUC_BT_enet1, AUC_BT_ridge1, result_BT_adapt1)
  save(result_BT1, file="result_BT1.RData")
  result_BT2 = cbind(AUC_BT_lasso2, AUC_BT_enet2, AUC_BT_ridge2, result_BT_adapt2)
  save(result_BT2, file = "result_BT2.RData")
  result_BT3 = cbind(AUC_BT_lasso3, AUC_BT_enet3, AUC_BT_ridge3, result_BT_adapt3)
  save(result_BT3, file="result_BT3.RData")
  