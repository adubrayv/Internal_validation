## cross-validation ##
## source
source("lambda1setuning.R")

## Settings of cross-validation: K = 10,  repetitions = 10, alpha = 0.95/0.5/0.05
### 1) Repeated cross validation using lambda.1se tuned in model selection for each scenario and penalized method
### 2) Extraction of the mean of AUC(t), Q10, Q90

##scenario1##
  ##lasso-like
  CV_lasso1 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data100[[i]], 
             0.95, 
             bestlambda1[[i]], 
             15000, 
             10,
             Weight = FALSE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_lasso1 = Get_AUC_CV(CV_lasso1)
  save(AUC_CV_lasso1, file = "result_CV_lasso1.RData")
  
  ##enet
  CV_enet1 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data100[[i]], 
             0.5, 
             bestlambdaenet1[[i]], 
             15000, 
             10,
             Weight = FALSE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_enet1 = Get_AUC_CV(CV_enet1)
  save(AUC_CV_enet1, file = "result_CV_enet1.RData")
  
  ##ridge-like
  CV_ridge1 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data100[[i]], 
             0.05, 
             bestlambdaridge1[[i]], 
             15000, 
             10,
             Weight = FALSE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_ridge1 = Get_AUC_CV(CV_ridge1)
  save(AUC_CV_ridge1, file = "result_CV_ridge1.RData")
  
  ##Adaptative Enet
  AUC_CV_adapt1 = lapply(1:100, function (i) CV_AUC(sim_data100[[i]], 
                                                      0.5, 
                                                      bestlambdaenet1[[i]], 
                                                      15000, 
                                                      10, 
                                                      Weight = TRUE))
  result_CV_adapt1 = Get_AUC_CV(AUC_CV_adapt1)
  save(result_CV_adapt1, file = "result_CV_adapt1.RData")

    ##scenario2##
    ##lasso-like
  CV_lasso2 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data500[[i]], 
             0.95, 
             bestlambda2[[i]], 
             15000, 
             10,
             Weight = TRUE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_lasso2 = Get_AUC_CV(CV_lasso2)
  save(AUC_CV_lasso2, file = "result_CV_lasso2.RData")
  
  ##enet
  CV_enet2 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data500[[i]], 
             0.5, 
             bestlambdaenet2[[i]], 
             15000, 
             10,
             Weight = TRUE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  
  AUC_CV_enet2 = Get_AUC_CV(CV_enet2)
  save(AUC_CV_enet2, file = "result_CV_enet2.RData")
  
  ##ridge-like
  CV_ridge2 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data500[[i]], 
             0.05, 
             bestlambdaridge2[[i]], 
             15000, 
             10,
             Weight = TRUE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_ridge2 = Get_AUC_CV(CV_ridge2)
  save(AUC_CV_ridge2, file = "result_CV_ridge2.RData")
  
  #adaptative Enet
  AUC_CV_adapt2 = lapply(1:100, function (i) CV_AUC(sim_data500[[i]], 
                                                    0.5, 
                                                    bestlambdaenet2[[i]],
                                                    15000, 
                                                    10, 
                                                    Weight = TRUE))
  result_CV_adapt2 = Get_AUC_CV(AUC_CV_adapt2)
  save(result_CV_adapt2, file = "result_CV_adapt2.RData")

##scenario3##
  ##lasso-like
  CV_lasso3 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data1000[[i]], 
             0.95, 
             bestlambda3[[i]], 
             15000, 
             10,
             Weight = TRUE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_lasso3 = Get_AUC_CV(CV_lasso3)
  save(AUC_CV_lasso3, file = "result_CV_lasso3.RData")
  
  ##enet
  CV_enet3 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data1000[[i]], 
             0.5, 
             bestlambdaenet3[[i]], 
             15000, 
             10,
             Weight = TRUE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_enet3 = Get_AUC_CV(CV_enet3)
  save(AUC_CV_enet3, file = "result_CV_enet3.RData")
  
  ##ridge-like
  CV_ridge3 <- lapply(1:100, function(i) {
    tryCatch({
      CV_AUC(sim_data1000[[i]], 
             0.05, 
             bestlambdaridge3[[i]], 
             15000, 
             10,
             Weight = TRUE)
    }, error = function(e) {
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # Return NULL or some default value
    })
  })
  AUC_CV_ridge3 = Get_AUC_CV(CV_ridge3)
  save(AUC_CV_ridge3, file = "result_CV_ridge3.RData")
  
  ##Adaptative Enet
  AUC_CV_adapt3 = lapply(1:100, function (i) CV_AUC(sim_data1000[[i]], 
                                                    0.5, 
                                                    bestlambdaenet3[[i]], 
                                                    15000, 
                                                    10, 
                                                    Weight = TRUE))
  result_CV_adapt3 = Get_AUC_CV(AUC_CV_adapt3)
  save(result_CV_adapt3, file = "result_CV_adapt3.RData")

    ## Save the results
  result_CV1 = cbind(AUC_CV_lasso1, AUC_CV_enet1, AUC_CV_ridge1, result_CV_adapt1)
  save(result_CV1, file = "result_CV1.RData")
  result_CV2 = cbind(AUC_CV_lasso2, AUC_CV_enet2, AUC_CV_ridge2, result_CV_adapt2)
  save(result_CV2, file = "result_CV2.RData")
  result_CV3 = cbind(AUC_CV_lasso3, AUC_CV_enet3, AUC_CV_ridge3, result_CV_adapt3)
  save(result_CV3, file = "result_CV3.RData")
  