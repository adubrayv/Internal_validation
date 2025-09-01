## Settings of cross-validation: K = 10,  repetitions = 10, alpha = 0.95/0.5/0.05
### 1) Repeated nested cross validation using lambda.1se tuned in model selection for each scenario and penalized regularization
### 2) Extraction of the mean of AUC(t)

NCV_AUC <- function(v, a, lambdas, n_genes, n_repeats = 5, Weight = TRUE, seed = 1234) {
  # Load necessary libraries
  library(survival)
  library(glmnet)
  library(timeROC)
  library(pec)
  library(survcomp)  # for concordance.index
  
  set.seed(seed)
  n <- nrow(v)
  ind_cov <- 1:(n_genes + 4)
  
  # Prepare data
  X <- as.matrix(v[, ind_cov])
  Y <- Surv(v$surv, v$event)
  
  # Step 1: Conditionally fit initial Cox model for adaptive weights
  if (Weight) {
    initial_cox_model <- glmnet(X, Y, alpha = 0.5, family = "cox", lambda = 0)
    initial_coefs <- as.numeric(coef(initial_cox_model))
    adaptive_weights <- 1 / (abs(initial_coefs) + 1e-4)
    penalty_factors <- adaptive_weights
  } else {
    penalty_factors <- rep(1, ncol(X))
  }
  
  # Store results across repetitions
  all_time_roc <- list()
  all_ibs <- list()
  all_cindex <- numeric()
  
  for (rep in 1:n_repeats) {
    # Create outer folds (5-fold CV)
    outer_folds <- createFolds(v$surv, k = 5, list = TRUE)
    
    # Prepare containers for each repetition - we only evaluate on the 3rd fold test here
    time_roc_rep <- list()
    ibs_rep <- list()
    
    # Hyperparameter tuning on first two outer folds combined
    tuning_indices <- unlist(outer_folds[1:4])
    test_fold_indices <- outer_folds[[5]]
    
    X_train_tuning <- X[tuning_indices, ]
    Y_train_tuning <- Y[tuning_indices]
    All_train_tuning <- v[tuning_indices, ]
    
    # Inner CV folds for tuning
    inner_folds <- createFolds(All_train_tuning$surv, k = 5, list = TRUE)
    cv_errors <- numeric(length(lambdas))
    
    for (lambda_idx in seq_along(lambdas)) {
      lambda <- lambdas[lambda_idx]
      inner_time_roc <- list()
      
      for (j in seq_along(inner_folds)) {
        inner_test_indices <- inner_folds[[j]]
        inner_train_indices <- setdiff(1:length(Y_train_tuning), inner_test_indices)
        
        X_inner_train <- X_train_tuning[inner_train_indices, ]
        Y_inner_train <- Y_train_tuning[inner_train_indices]
        All_inner_train <- All_train_tuning[inner_train_indices, ]
        
        X_inner_val <- X_train_tuning[inner_test_indices, ]
        Y_inner_val <- Y_train_tuning[inner_test_indices]
        All_inner_val <- All_train_tuning[inner_test_indices, ]
        
        model_inner <- glmnet(X_inner_train, Y_inner_train, alpha = a, family = "cox",
                              lambda = lambda, penalty.factor = penalty_factors)
        
        All_inner_val$preds <- predict(model_inner, newx = X_inner_val, s = lambda, type = "link")
        
        inner_time_roc[[j]] <- timeROC(T = All_inner_val$surv, delta = All_inner_val$event, marker = All_inner_val$preds,
                                       weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
      }
      
      cv_errors[lambda_idx] <- mean(unlist(inner_time_roc), na.rm = TRUE)
    }
    
    bestlambda <- lambdas[which.max(cv_errors)]
    
    # Fit model on combined tuning folds using best lambda
    model <- glmnet(X_train_tuning, Y_train_tuning, alpha = a, family = "cox",
                    lambda = bestlambda, penalty.factor = penalty_factors)
    
    # Outer fold test data
    X_test <- X[test_fold_indices, ]
    Y_test <- Y[test_fold_indices]
    All_test <- v[test_fold_indices, ]
    
    All_test$preds <- predict(model, newx = X_test, s = bestlambda, type = "link")
    
    if (any(is.na(All_test$preds))) {
      cox_formula <- as.formula(Surv(surv, event) ~ 1)
      marker <- 1
      c_index <- NA
    } else {
      cox_formula <- as.formula(Surv(surv, event) ~ preds)
      marker <- All_test$preds
      
      # Calculate c-index using survcomp
      c_index <- concordance.index(
        x = All_test$preds,
        surv.time = All_test$surv,
        surv.event = All_test$event
      )$c.index
      
      # Alternatively, using survival package only:
      # c_index <- survConcordance(Surv(All_test$surv, All_test$event) ~ All_test$preds)$concordance
    }
    
    # Calculate time-dependent AUC
    time_roc_rep[[1]] <- timeROC(T = All_test$surv, delta = All_test$event, marker = marker,
                                 weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
    
    # Fit Cox model for IBS calculation
    cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
    
    if (any(is.na(coef(cox_fit)))) {
      next  # Skip iteration if NA coefficients
    }
    
    pecfit <- pec(object = list(cox_fit), formula = Surv(surv, event) ~ preds, data = All_test,
                  exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
    
    ibs_rep[[1]] <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0)
    
    # Store results
    all_time_roc[[rep]] <- time_roc_rep[[1]]
    all_ibs[[rep]] <- ibs_rep[[1]]
    all_cindex[rep] <- c_index
  }
  
  # Summary across repetitions
  time_roc <- colMeans(as.data.frame(do.call(rbind, all_time_roc)), na.rm = TRUE)
  ibs_mean <- colMeans(as.data.frame(do.call(rbind, all_ibs)), na.rm = TRUE)
  c_index_mean <- mean(all_cindex, na.rm = TRUE)
  
  return(list(time_roc = time_roc, ibs = ibs_mean, c_index = c_index_mean, bestmodel = model))
}

##as examples or lambda.grid to define
lambda_lasso = seq(0,0.35,0.01)
lambda_enet = seq(0,0.7,0.01)
lambda_ridge = seq(0.3,6,0.1)
