## cross-validation ## using apparent approach, development of the model with fixed hyperparameters before the assessment of predictive performances
## source
source("lambda1setuning.R")

## Settings of cross-validation: K = 10,  repetitions = 10, alpha = 0.95/0.5/0.05
### 1) Repeated cross validation using lambda.1se tuned in model selection for each scenario and penalized method
### 2) Extraction of the mean of AUC(t)
CV_AUC <- function(v, a, bestlambda, n_genes, n_repeats = 10, Weight = TRUE, seed = 1234) {
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
    adaptive_weights <- 1 / (abs(initial_coefs) + 1e-4) 
    penalty_factors <- adaptive_weights
  } else {
    penalty_factors <- rep(1, ncol(X))
  }
  
  # Store results across repetitions
  all_time_roc <- list()
  all_ibs <- list()
  all_cindex <- list()
  
  for (rep in 1:n_repeats) {
    folds <- createFolds(v$surv, k = 10, list = TRUE)
    time_roc_rep <- list()
    ibs_rep <- list()
    cindex_rep <- numeric(length(folds))
    
    for (i in seq_along(folds)) {
      test_indices <- folds[[i]]
      train_indices <- setdiff(1:n, test_indices)
      
      # Training data
      X_train <- X[train_indices, ]
      Y_train <- Y[train_indices]
      All_train <- v[train_indices, ]
      
      # Test data
      X_test <- X[test_indices, ]
      All_test <- v[test_indices, ]
      
      # Fit the model with weighted elastic net (adaptive penalty)
      fit_model <- tryCatch({
        glmnet(X_train, Y_train, alpha = a, family = "cox", lambda = bestlambda, penalty.factor = penalty_factors)
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(fit_model)) {
        next  # Skip this fold if fitting failed
      }
      
      # Make predictions on the test data
      All_test$preds <- predict(fit_model, newx = X_test, s = bestlambda, type = "link")
      if (any(is.na(All_test$preds))) {
        next  # Skip fold if NA in predictions
      }
      
      # Calculate time-dependent AUC for the fold
      time_roc_rep[[i]] <- timeROC(T = All_test$surv, delta = All_test$event, marker = All_test$preds, 
                                   weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
      
      # Fit Cox model for IBS and c-index calculation
      cox_formula <- as.formula(Surv(surv, event) ~ preds)
      cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
      
      if (any(is.na(coef(cox_fit)))) {
        ibs_rep[[i]] <- rep(NA, length(seq(0, 4.9, 0.1)))
        cindex_rep[i] <- NA
        next
      }
      
      # Compute IBS using PEC
      pecfit <- pec(object = list(cox_fit), formula = Surv(surv, event) ~ preds, data = All_test,
                    exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
      ibs_rep[[i]] <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0)
      
      # Calculate c-index using survcomp package
      cindex_rep[i] <- concordance.index(
        x = All_test$preds,
        surv.time = All_test$surv,
        surv.event = All_test$event
      )$c.index
    }
    
    # Aggregate results per repetition
    all_time_roc[[rep]] <- colMeans(as.data.frame(do.call(rbind, time_roc_rep)), na.rm = TRUE)
    all_ibs[[rep]] <- colMeans(as.data.frame(do.call(rbind, ibs_rep)), na.rm = TRUE)
    all_cindex[[rep]] <- mean(cindex_rep, na.rm = TRUE)  # Mean c-index across folds
  }
  
  return(list(
    time_roc = all_time_roc,
    ibs = all_ibs,
    c_index = all_cindex
  ))
}

