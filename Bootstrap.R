## source
source("lambda1setuning.R")
##library
library(parallel)
library(foreach)
library(doParallel)
library(future)
library(furrr)
## Settings of bootstrap: original B = 100, a = 0.95/0.5/0.05
### 1) Original bootstrap using lambda.1se tuned in model selection for each scenario and penalized method
### approach of conventional boostrap used in HDNOM package without correction
### 2) Extraction of the mean of AUC(t)
AUC_bootstrap <- function(v, a, bestlambda, n_genes, n_bootstraps = 100, Weight = TRUE, seed = 1234) {
  # Load necessary libraries
  library(survival)
  library(glmnet)
  library(timeROC)
  library(pec)
  library(survcomp)
  
  set.seed(seed)
  n <- nrow(v)
  ind_cov <- 1:(n_genes + 4)
  
  # Prepare data
  X <- as.matrix(v[, ind_cov])
  Y <- Surv(v$surv, v$event)
  
  # Step 1: Conditionally fit an initial Cox model for adaptive weights if Weight = TRUE
  if (Weight) {
    initial_cox_model <- glmnet(X, Y, alpha = 0.5, family = "cox", lambda = 0)
    initial_coefs <- as.numeric(coef(initial_cox_model))  # Extract the initial coefficients
    # Step 2: Compute adaptive weights (inverse of absolute coefficient values)
    adaptive_weights <- 1 / (abs(initial_coefs) + 1e-4)
    # Step 3: Set penalty factors accordingly
    penalty_factors <- adaptive_weights
  } else {
    penalty_factors <- rep(1, ncol(X))
  }
  
  # Store results
  time_roc <- list()
  ibs <- list()
  c_indices <- numeric()  # For c-index across bootstraps
  
  # Bootstrap iterations
  for (i in 1:n_bootstraps) {
    # Generate bootstrap sample indices
    train_indices <- sample(1:n, size = n, replace = TRUE)
    test_indices <- setdiff(1:n, train_indices)
    
    # Ensure there are out-of-bag (OOB) samples for testing
    if (length(test_indices) == 0) next
    
    # Training data
    X_train <- X[train_indices, ]
    Y_train <- Y[train_indices]
    All_train <- v[train_indices, ]
    
    # Test data (for prediction and evaluation)
    X_test <- X
    Y_test <- Y
    All_test <- v
    
    # Fit model with adaptive elastic net
    model <- glmnet(X_train, Y_train, alpha = a, family = "cox", lambda = bestlambda, penalty.factor = penalty_factors)
    
    # Make predictions on the test data
    All_test$preds <- predict(model, newx = X_test, s = bestlambda, type = "link")
    
    if (any(is.na(All_test$preds))) {
      cox_formula <- as.formula(Surv(surv, event) ~ 1)
      marker <- 1
      c_index <- NA
    } else {
      cox_formula <- as.formula(Surv(surv, event) ~ preds)
      marker <- All_test$preds
      
      # Calculate c-index (using survcomp)
      c_index <- concordance.index(
        x = All_test$preds,
        surv.time = All_test$surv,
        surv.event = All_test$event
      )$c.index
    }
    
    # Calculate time-dependent AUC
    time_roc[[i]] <- timeROC(T = All_test$surv, delta = All_test$event, marker = marker, 
                             weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
    
    # Fit Cox model for IBS calculation
    cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
    
    # Check for NA coefficients and handle
    if (any(is.na(coef(cox_fit)))) {
      return(list(time_roc = NA, ibs = NA, c_index = NA))
    }
    
    # Calculate Integrated Brier Score (IBS)
    pecfit <- pec(object = list(cox_fit), formula = cox_formula, data = All_test,
                  exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
    ibs[[i]] <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0)
    
    # Store c-index for this bootstrap
    c_indices[i] <- c_index
  }
  
  # Combine bootstrap results (averages)
  time_roc_mean <- colMeans(as.data.frame(do.call(rbind, time_roc)), na.rm = TRUE)
  ibs_mean <- colMeans(as.data.frame(do.call(rbind, ibs)), na.rm = TRUE)
  c_index_mean <- mean(c_indices, na.rm = TRUE)
  
  return(list(time_roc = time_roc_mean, ibs = ibs_mean, c_index = c_index_mean))
}
  

