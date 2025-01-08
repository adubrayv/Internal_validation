##Nested cross-validation ##
library(parallel)
library(foreach)
library(doParallel)
library(fastmap)
## Settings of cross-validation: K = 5,  repetitions = 5, alpha = 0.95/0.5/0.05
### 1) Repeated cross validation using lambda.1se tuned in model selection for each scenario and penalized method
### 2) Extraction of the mean of AUC(t), Q10, Q90
lambda_lasso = seq(0,0.35,0.01)
lambda_enet = seq(0,0.7,0.01)
lambda_ridge = seq(0.3,6,0.1)

##scenario1##
  ##lasso-like
  NCV_lasso1 = lapply(seq_along(sim_data100), function (i) NCV_AUC(sim_data100[[i]],
                                                                  0.95,
                                                                  lambdas = lambda_lasso,
                                                                  15000,
                                                                  5))
  
  AUC_NCV_lasso1 = Get_AUC_NCV(NCV_lasso1)

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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data100", "NCV_AUC"))
  NCV_enet1 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data100[[i]],
            0.5,
            lambdas = seq(0,0.7,0.01),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_enet1 = Get_AUC_NCV(NCV_enet1)
  
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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data100", "NCV_AUC"))
  NCV_enet1 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data100[[i]],
            0.05,
            lambdas = seq(0.3,6,0.1),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_ridge1 = Get_AUC_NCV(NCV_enet1)
  
##scenario2##
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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data500", "NCV_AUC"))
  NCV_lasso2 <- parLapply(cl, 1:100, function(i) {
     NCV_AUC(sim_data500[[i]],
             0.95,
             lambdas = seq(0,0.35,0.01),
             15000,
             5)
  })
  stopCluster(cl)
  AUC_NCV_lasso2 = Get_AUC_NCV(NCV_lasso2)

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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data500", "NCV_AUC"))
  NCV_enet2 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data500[[i]],
            0.3,
            lambdas = seq(0,0.7,0.01),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_enet2 = Get_AUC_NCV(NCV_enet2)
  
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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data500", "NCV_AUC"))
  NCV_ridge2 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data500[[i]],
            0.05,
            lambdas = seq(0.3,6,0.1),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_ridge2 = Get_AUC_NCV(NCV_ridge2)
  
##scenario3##
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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data1000", "NCV_AUC"))
  NCV_lassobis3 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data1000[[i]],
            0.95,
            lambdas = seq(0.02,0.20,0.01),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_lasso3 = Get_AUC_NCV(NCV_lasso3)
  
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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data1000", "NCV_AUC"))
  NCV_enet3 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data1000[[i]],
            0.5,
            lambdas = seq(0,0.35,0.01),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_enet3 = Get_AUC_NCV(NCV_enet3)
  
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
    library(caret)
  })
  # Export necessary objects and functions
  parallel::clusterExport(cl, varlist = c("sim_data1000", "NCV_AUC"))
  NCV_ridge3 <- parLapply(cl, 1:100, function(i) {
    NCV_AUC(sim_data1000[[i]],
            0.05,
            lambdas = seq(0.3,4,0.1),
            15000,
            5)
  })
  stopCluster(cl)
  AUC_NCV_ridge3 = Get_AUC_NCV(NCV_ridge3)
  
    ## Save the results
  result_NCV1 = cbind(AUC_NCV_lasso1, AUC_NCV_enet1, AUC_NCV_ridge1)
  save(result_NCV1, file = "result_NCV1.RData")
  result_NCV2 = cbind(AUC_NCV_lasso2, AUC_NCV_enet2, AUC_NCV_ridge2)
  save(result_NCV2, file = "result_NCV2.RData")
  result_NCV3 = cbind(AUC_NCV_lasso3, AUC_NCV_enet3, AUC_NCV_ridge3)
  save(result_NCV3, file = "result_NCV3.RData")
  
