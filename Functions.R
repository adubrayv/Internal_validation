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
                    
##split_sample validation

## sample division of datasets using 2/3 of training sample and 1/3 of testing sample
###prediction with mean AUC(t), Q10 and Q90 on test sample and integrated Brier Score (t)

split_sample_AUC <- function(v, a, bestlambda, n_genes, Weight = TRUE) {
  set.seed(1234)
  n <- nrow(v)
  # Split the data into training and testing sets
  ind_train <- sample(1:n, size = floor(0.667 * n))
  ind_test <- setdiff(1:n, ind_train)
  ind_cov <- 1:(n_genes + 4)  # Covariates include gene variables and 4 additional ones
  # Define the training set
  X_train <- as.matrix(v[ind_train, ind_cov])
  Y_train <- Surv(v$surv[ind_train], v$event[ind_train])
  # Define the testing set
  X_test <- as.matrix(v[ind_test, ind_cov])
  All_test <- v[ind_test, ]
  # Step 1: Conditionally fit an initial Cox model for adaptive weights if Weight = TRUE
  if (Weight) {
    initial_cox_model <- glmnet(X_train, Y_train, alpha = 0, family = "cox", lambda = 0)
    initial_coefs <- as.numeric(coef(initial_cox_model))  # Extract the initial coefficients
    # Step 2: Compute adaptive weights (inverse of absolute coefficient values)
    adaptive_weights <- 1 / (abs(initial_coefs) + 1e-4)  # Add small constant to avoid division by zero
    # Step 3: Set penalty factors to 0 for variables 15001 to 15004 (indices n_genes + 1 to n_genes + 4)
    penalty_factors <- adaptive_weights
    penalty_factors[(n_genes + 1):(n_genes + 4)] <- 0  # Unpenalized variables
  } else {
    # If Weight = FALSE, all variables have equal penalty
    penalty_factors <- rep(1, ncol(X_train))
  }
  # Step 4: Fit the Cox model with adaptive elastic net using the computed penalty factors
  model <- glmnet(X_train, Y_train, alpha = a, family = "cox", lambda = bestlambda, penalty.factor = penalty_factors)
  # Predict on the testing set
  All_test$preds <- predict(model, newx = X_test, s = bestlambda, type = "link")
  if (any(is.na(All_test$preds))) {
    # Use an intercept-only model if any NA is present in predictions
    cox_formula <- as.formula(Surv(surv, event) ~ 1)
    marker <- 1
  } else {
    # Use the full model if no NA is present in predictions
    cox_formula <- as.formula(Surv(surv, event) ~ preds)
    marker <- All_test$preds
  }
  # Calculate AUC using time-dependent ROC
  AUCroc <- timeROC(T = All_test$surv, delta = All_test$event, marker = marker, cause = 1,
                    times = seq(0.9, 4.9, 0.1), ROC = TRUE)$AUC
  # Fit Cox proportional hazards model
  cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
  if (any(is.na(coef(cox_fit)))) {
    # Return NA if coefficients contain NA
    return(list(time_roc = as.numeric(rep(NA, length(seq(0.9, 4.9, 0.1)))),
                ibs = as.numeric(rep(NA, length(seq(0, 4.9, 0.1))))))
  }
  # Calculate Integrated Brier Score using PEC
  pecfit <- pec(object = list(cox_fit), formula = Surv(surv, event) ~ preds, data = All_test,
                exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
  # Store the Integrated Brier Score (IBS)
  ibs <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0.1)
  return(list(time_roc = unname(AUCroc), ibs = unname(ibs[1, 1:50])))
}

## Extraction of AUC(t) and integrated IBS(t) of split_sample process with mean, Q10 and Q90 values
Get_AUC_SSV <- function(split_sample_AUC) {
  # Initialize lists to store results from all folds
  all_time_roc <- list()
  all_ibs_mean <- list()
  
  for (i in seq_along(split_sample_AUC)) {
    # Only include results without NA values
    all_time_roc[[i]] <- as.data.frame(split_sample_AUC[[i]]$time_roc)
    all_ibs_mean[[i]] <- as.data.frame(split_sample_AUC[[i]]$ibs)
  }
  AUC_SSV <- as.data.frame(do.call(cbind, all_time_roc))
  ibs_SSV <- as.data.frame(do.call(cbind, all_ibs_mean))
  # Calculate summary statistics for time_roc
  time_roc_df <- data.frame(
    Time = seq(0.9, 4.9, 0.1),  # Assuming these are the time points used
    Mean = apply(AUC_SSV,1, mean, na.rm = TRUE),
    Q10 = apply(AUC_SSV, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(AUC_SSV, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  # Calculate summary statistics for ibs_mean
  ibs_mean_df <- data.frame(
    Time = seq(0, 4.9, 0.1),
    Mean = apply(ibs_SSV,1, mean, na.rm = TRUE),
    Q10 = apply(ibs_SSV, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(ibs_SSV, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  return(list(time_roc_summary = time_roc_df, ibs_mean_summary = ibs_mean_df))
}

## Bootstrap conventional
AUC_bootstrap <- function(v, a, bestlambda, n_genes, n_bootstraps = 100, Weight = TRUE) {
  set.seed(1234)
  n <- nrow(v)
  ind_cov <- 1:(n_genes + 4)
  # Prepare data
  X <- as.matrix(v[, ind_cov])
  Y <- Surv(v$surv, v$event)
  # Step 1: Conditionally fit an initial Cox model for adaptive weights if Weight = TRUE
  if (Weight) {
    initial_cox_model <- glmnet(X, Y, alpha = 0, family = "cox", lambda = 0)
    initial_coefs <- as.numeric(coef(initial_cox_model))  # Extract the initial coefficients
    # Step 2: Compute adaptive weights (inverse of absolute coefficient values)
    adaptive_weights <- 1 / (abs(initial_coefs) + 1e-4)  # Add small constant to avoid division by zero
    # Step 3: Set penalty factors to 0 for variables 15001 to 15004 (indices n_genes + 1 to n_genes + 4)
    penalty_factors <- adaptive_weights
    penalty_factors[(n_genes + 1):(n_genes + 4)] <- 0  # Unpenalized variables
  } else {
    # If Weight = FALSE, all variables have equal penalty
    penalty_factors <- rep(1, ncol(X))
  }
  # Store results
  time_roc <- list()
  pecfit <- list()
  ibs <- list()
  # Bootstrap iterations
  for (i in 1:n_bootstraps) {
    # Generate bootstrap sample indices
    train_indices <- sample(1:n, size = n, replace = TRUE)
    test_indices <- setdiff(1:n, train_indices)
    # Ensure there are out-of-bag (OOB) samples
    if (length(test_indices) == 0) next
    # Training data
    X_train <- X[train_indices, ]
    Y_train <- Y[train_indices]
    All_train <- v[train_indices, ]
    # Test data (using all data for prediction in bootstrap)
    X_test <- X
    Y_test <- Y
    All_test <- v
    # Fit the model with weighted elastic net (adaptive penalty and unpenalized variables)
    model <- glmnet(X_train, Y_train, alpha = a, family = "cox", lambda = bestlambda, penalty.factor = penalty_factors)
    # Make predictions on the test data
    All_test$preds <- predict(model, newx = X_test, s = bestlambda, type = "link")
    # Check for NA values in predictions
    if (any(is.na(All_test$preds))) {
      # Replace ~ preds with ~1 if any NA is present
      cox_formula <- as.formula(Surv(surv, event) ~ 1)
      marker <- 1
    } else {
      # Use ~preds if no NA is present
      cox_formula <- as.formula(Surv(surv, event) ~ preds)
      marker <- All_test$preds
    }
    # Calculate AUC for the predictions of this fold
    time_roc[[i]] <- timeROC(T = All_test$surv, delta = All_test$event, marker = marker, 
                             weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
    # Fit a Cox model using the chosen formula
    cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
    # Check for NA in coefficients and handle it
    if (any(is.na(coef(cox_fit)))) {
      return(list(time_roc = NA, ibs = NA))
    }
    # Calculate the Integrated Brier Score using PEC
    pecfit <- pec(object = list(cox_fit), formula = cox_formula, data = All_test,
                  exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
    
    # Store the Integrated Brier Score (IBS)
    ibs[[i]] <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0)
  }
  # Combine all AUC and IBS results into single data frames
  time_roc <- colMeans(as.data.frame(do.call(rbind, time_roc)), na.rm = TRUE)
  ibs_mean <- colMeans(as.data.frame(do.call(rbind, ibs)), na.rm = TRUE)
  return(list(time_roc = time_roc, ibs = ibs_mean))
}

## Extraction of AUC(t) and integrated IBS(t) of conventional bootstrap process with Mean, Q10 and Q90 values
Get_AUC_BT <- function(AUC_bootstrap) {
  # Initialize lists to store results from all folds
  all_time_roc <- list()
  all_ibs_mean <- list()
  # Loop over each element of CVtest_AUC to aggregate time_roc and ibs_mean results
  for (i in seq_along(AUC_bootstrap)) {
    # Combine results into data frames
    all_time_roc[[i]] <- as.data.frame(AUC_bootstrap[[i]]$time_roc)
    all_ibs_mean[[i]] <- as.data.frame(AUC_bootstrap[[i]]$ibs)
  }
  # Combine all results into single data frames
  AUC_BT <- as.data.frame(do.call(cbind, all_time_roc))
  ibs_BT <- as.data.frame(do.call(cbind, all_ibs_mean))
  # Calculate summary statistics for time_roc
  time_roc_df <- data.frame(
    Time = seq(1, 4.9, 0.1),  # Assuming these are the time points used
    Mean = rowMeans(AUC_BT, na.rm = TRUE),
    Q10 = apply(AUC_BT, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(AUC_BT, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  # Calculate summary statistics for ibs_mean
  ibs_mean_df <- data.frame(
    Time = seq(0, 4.9, 0.1),
    Mean = rowMeans(ibs_BT, na.rm = TRUE),
    Q10 = apply(ibs_BT, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(ibs_BT, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  return(list(time_roc_summary = time_roc_df, ibs_mean_summary = ibs_mean_df))
}

##repeated k fold cross_validation

CV_AUC <- function(v, a, bestlambda, n_genes, n_repeats = 10) {
  set.seed(1234)
  n <- nrow(v)
  ind_cov <- 1:(n_genes + 4)
  # Prepare data
  X <- as.matrix(v[, ind_cov])
  Y <- Surv(v$surv, v$event)
  # Store results across repetitions
  all_time_roc <- list()
  all_ibs <- list()
  # Repeat the cross-validation process n_repeats times
  for (rep in 1:n_repeats) {
    # Create folds for cross-validation
    folds <- createFolds(v$surv, k = 10, list = TRUE)
    # Store results for each fold within the current repetition
    time_roc_rep <- list()
    ibs_rep <- list()
    # Loop over each fold
    for (i in seq_along(folds)) {
      test_indices <- folds[[i]]
      train_indices <- setdiff(1:n, test_indices)
      # Training data
      X_train <- X[train_indices, ]
      Y_train <- Y[train_indices]
      All_train <- v[train_indices, ]
      # Test data
      X_test <- X[test_indices, ]
      Y_test <- Y[test_indices]
      All_test <- v[test_indices, ]
      fit_model <- tryCatch({
        glmnet(X_train, Y_train, alpha = a, family = "cox", lambda = bestlambda)
      }, error = function(e) {
        return(NULL)  # Return NULL if model fitting fails
      })
      if (is.null(fit_model)) {
        next  # Skip this fold if model fitting failed
      }
      # Make predictions on the test data
      All_test$preds <- predict(fit_model, newx = X_test, s = bestlambda, type = "link")
      if (any(is.na(All_test$preds))) {
        next  # Skip to the next iteration of the fold loop
      }
      # Calculate AUC for the predictions of this fold
      time_roc_rep[[i]] <- timeROC(T = All_test$surv, delta = All_test$event, marker = All_test$preds, 
                                   weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
      # Fit a Cox model using the predicted linear predictors
      cox_formula <- as.formula(Surv(surv, event) ~ preds)
      cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
      if (any(is.na(coef(cox_fit)))) {
        ibs_rep[[i]] <- rep(NA, length(seq(0, 4.9, 0.1)))  # Store NA if coefficients contain NA
        next
      }
      # Calculate the Integrated Brier Score using PEC
      pecfit <- pec(object = list(cox_fit), formula = Surv(surv, event) ~ preds, data = All_test,
                    exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
      # Store the Integrated Brier Score (IBS)
      ibs_rep[[i]] <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0)
    }
    # Combine results from all folds within the current repetition
    all_time_roc[[rep]] <- colMeans(as.data.frame(do.call(rbind, time_roc_rep)), na.rm = TRUE)
    all_ibs[[rep]] <- colMeans(as.data.frame(do.call(rbind, ibs_rep)), na.rm = TRUE)
  }
  return(list(time_roc = all_time_roc, ibs = all_ibs))
}

Get_AUC_CV <- function(CV_AUC) {
  # Initialize lists to store results from all folds
  all_time_roc <- list()
  all_ibs_mean <- list()
  # Loop over each element of CV_AUC to aggregate time_roc and ibs_mean results
  for (i in seq_along(CV_AUC)) {
    # Check if the current element is NULL
    if (is.null(CV_AUC[[i]])) {
      # Add NA values for time_roc and ibs_mean if NULL
      all_time_roc[[i]] <- as.data.frame(rep(NA, length(seq(1, 4.9, 0.1))))
      all_ibs_mean[[i]] <- as.data.frame(rep(NA, length(seq(0, 4.9, 0.1))))
      next  # Skip to the next iteration
    }
    # Combine results into data frames
    all_time_roc[[i]] <- as.data.frame(CV_AUC[[i]]$time_roc)
    all_ibs_mean[[i]] <- as.data.frame(CV_AUC[[i]]$ibs)
  }
  # Combine all non-NULL results into single data frames
  AUC_CV <- as.data.frame(do.call(cbind, all_time_roc))
  ibs_CV <- as.data.frame(do.call(cbind, all_ibs_mean))
  # Calculate summary statistics for time_roc
  time_roc_df <- data.frame(
    Time = seq(1, 4.9, 0.1),  # Assuming these are the time points used
    Mean = rowMeans(AUC_CV, na.rm = TRUE),
    Q10 = apply(AUC_CV, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(AUC_CV, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  # Calculate summary statistics for ibs_mean
  ibs_mean_df <- data.frame(
    Time = seq(0, 4.9, 0.1),
    Mean = rowMeans(ibs_CV, na.rm = TRUE),
    Q10 = apply(ibs_CV, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(ibs_CV, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  return(list(time_roc_summary = time_roc_df, ibs_mean_summary = ibs_mean_df))
}


##nested Cross Validation defined by inner and outer loop
NCV_AUC <- function(v, a, lambdas, n_genes, n_repeats = 10, Weight = TRUE) {
  set.seed(1234)
  n <- nrow(v)
  ind_cov <- 1:(n_genes + 4)
  # Prepare data
  X <- as.matrix(v[, ind_cov])
  Y <- Surv(v$surv, v$event)
  # Store results across repetitions
  all_time_roc <- list()
  all_ibs <- list()
  # Repeat the cross-validation process n_repeats times
  for (rep in 1:n_repeats) {
    # Create folds for the outer cross-validation
    outer_folds <- createFolds(v$surv, k = 5, list = TRUE)
    # Store results for each repetition
    time_roc_rep <- list()
    ibs_rep <- list()
    # Perform hyperparameter tuning on the first four outer folds
    tuning_indices <- unlist(outer_folds[1:4])
    test_fold_indices <- outer_folds[[5]]
    # Prepare the tuning dataset (first four folds)
    X_train_tuning <- X[tuning_indices, ]
    Y_train_tuning <- Y[tuning_indices]
    All_train_tuning <- v[tuning_indices, ]
    # Inner cross-validation for hyperparameter tuning (5-fold CV)
    inner_folds <- createFolds(All_train_tuning$surv, k = 5, list = TRUE)
    cv_errors <- numeric(length(lambdas))
    for (lambda_idx in seq_along(lambdas)) {
      lambda <- lambdas[lambda_idx]
      inner_time_roc <- list()
      # Loop over each inner fold
      for (j in seq_along(inner_folds)) {
        inner_test_indices <- inner_folds[[j]]
        inner_train_indices <- setdiff(1:length(Y_train_tuning), inner_test_indices)
        # Inner training and validation data
        X_inner_train <- X_train_tuning[inner_train_indices, ]
        Y_inner_train <- Y_train_tuning[inner_train_indices]
        All_inner_train <- All_train_tuning[inner_train_indices, ]
        X_inner_val <- X_train_tuning[inner_test_indices, ]
        Y_inner_val <- Y_train_tuning[inner_test_indices]
        All_inner_val <- All_train_tuning[inner_test_indices, ]
        # Fit the model with the current lambda using penalty factors
        model_inner <- glmnet(X_inner_train, Y_inner_train, alpha = a, family = "cox", lambda = lambda)
        # Predictions for the inner validation set
        All_inner_val$preds <- predict(model_inner, newx = X_inner_val, s = lambda, type = "link")
        # Compute time-dependent AUC for this inner fold
        inner_time_roc[[j]] <- timeROC(T = All_inner_val$surv, delta = All_inner_val$event, marker = All_inner_val$preds, 
                                       weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
      }
      # Compute mean error across inner folds for current lambda
      cv_errors[lambda_idx] <- mean(unlist(inner_time_roc), na.rm = TRUE)
    }
    # Select best lambda (highest AUC)
    bestlambda <- lambdas[which.max(cv_errors)]
    # Fit the model on the first four outer folds using the selected best lambda
    model <- glmnet(X_train_tuning, Y_train_tuning, alpha = a, family = "cox", lambda = bestlambda)
    # Test on the 5th outer fold
    X_test <- X[test_fold_indices, ]
    Y_test <- Y[test_fold_indices]
    All_test <- v[test_fold_indices, ]
    # Predictions for the outer test set
    All_test$preds <- predict(model, newx = X_test, s = bestlambda, type = "link")
    # Check for NA values in predictions
    if (any(is.na(All_test$preds))) {
      # Replace ~ preds with ~1 if any NA is present
      cox_formula <- as.formula(Surv(surv, event) ~ 1)
      marker <- 1 
    } else {
      # Use ~preds if no NA is present
      cox_formula <- as.formula(Surv(surv, event) ~ preds)
      marker <- All_test$preds
    }
    # Calculate AUC for the predictions of the 5th outer fold
    time_roc_rep[[1]] <- timeROC(T = All_test$surv, delta = All_test$event, marker = marker, 
                                 weighting = "marginal", cause = 1, times = seq(1, 4.9, 0.1), ROC = TRUE)$AUC
    # Fit a Cox model using the predicted linear predictors
    cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
    if (any(is.na(coef(cox_fit)))) {
      # Skip this iteration if NA coefficients are found
      next
    }
    # Calculate the Integrated Brier Score using PEC
    pecfit <- pec(object = list(cox_fit), formula = Surv(surv, event) ~ preds, data = All_test,
                  exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
    # Store the Integrated Brier Score (IBS)
    ibs_rep[[1]] <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0)
    # Combine results from the 5th outer fold
    all_time_roc[[rep]] <- time_roc_rep[[1]]
    all_ibs[[rep]] <- ibs_rep[[1]]
  }
  # Combine results across all repetitions using mean of last fold repetitions
  time_roc <- colMeans(as.data.frame(do.call(rbind, all_time_roc)), na.rm = TRUE)
  ibs_mean <- colMeans(as.data.frame(do.call(rbind, all_ibs)), na.rm = TRUE)
  return(list(time_roc = time_roc, ibs = ibs_mean, bestmodel = model))
}

## Extraction of nested_cross validation data: mean AUC and IBS obtained on last fold using lambda.1se hypertuning  on the four first inner loop
Get_AUC_NCV <- function(NCV_AUC) {
  # Initialize lists to store results from all folds
  all_time_roc <- list()
  all_ibs_mean <- list()
  # Loop over each element of NCV_AUC to aggregate time_roc and ibs_mean results
  for (i in seq_along(NCV_AUC)) {
    # Check if the current element or its time_roc or ibs fields are NULL or zero-length
    if (is.null(NCV_AUC[[i]]) || length(NCV_AUC[[i]]$time_roc) == 0 || length(NCV_AUC[[i]]$ibs) == 0) {
      # Add NA values for time_roc and ibs_mean if NULL or zero-length
      all_time_roc[[i]] <- as.data.frame(rep(NA, length(seq(1, 4.9, 0.1))))
      all_ibs_mean[[i]] <- as.data.frame(rep(NA, length(seq(0, 4.9, 0.1))))
      next  # Skip to the next iteration
    }
    # Combine results into data frames
    all_time_roc[[i]] <- as.data.frame(NCV_AUC[[i]]$time_roc)
    all_ibs_mean[[i]] <- as.data.frame(NCV_AUC[[i]]$ibs)
  }
  # Combine all results into single data frames
  AUC_NCV <- as.data.frame(do.call(cbind, all_time_roc))
  ibs_NCV <- as.data.frame(do.call(cbind, all_ibs_mean))
  # Calculate summary statistics for time_roc
  time_roc_df <- data.frame(
    Time = seq(1, 4.9, 0.1),  # Assuming these are the time points used
    Mean = rowMeans(AUC_NCV, na.rm = TRUE),
    Q10 = apply(AUC_NCV, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(AUC_NCV, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  # Calculate summary statistics for ibs_mean
  ibs_mean_df <- data.frame(
    Time = seq(0, 4.9, 0.1),
    Mean = rowMeans(ibs_NCV, na.rm = TRUE),
    Q10 = apply(ibs_NCV, 1, function(x) quantile(x, probs = 0.10, na.rm = TRUE)),
    Q90 = apply(ibs_NCV, 1, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  )
  return(list(time_roc_summary = time_roc_df, ibs_mean_summary = ibs_mean_df))
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



