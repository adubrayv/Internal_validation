##split sample
## source
source("lambda1setuning.R")

## sample division of datasets using 70% of training sample and 30% of testing sample
###optimization of lambda.1se on training sample and prediction with mean AUC(t)
split_sample_AUC <- function(v, a, bestlambda, n_genes, Weight = TRUE, seed = 1234) {
  # You may need to install 'survcomp' if not already done:
  # install.packages("survcomp")
  library(survival)
  library(glmnet)
  library(timeROC)
  library(pec)
  library(survcomp)
  
  set.seed(seed)
  n <- nrow(v)
  # Split the data into training and testing sets
  ind_train <- sample(1:n, size = floor(0.7 * n))
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
    initial_cox_model <- glmnet(X_train, Y_train, alpha = 0.5, family = "cox", lambda = 0)
    initial_coefs <- as.numeric(coef(initial_cox_model))  # Extract the initial coefficients
    # Step 2: Compute adaptive weights (inverse of absolute coefficient values)
    adaptive_weights <- 1 / (abs(initial_coefs) + 1e-4) 
    # Step 3: Set penalty factors to 0 for variables 15001 to 15004 (indices n_genes + 1 to n_genes + 4)
    penalty_factors <- adaptive_weights
  } else {
    # If Weight = FALSE, all variables have equal penalty
    penalty_factors <- rep(1, ncol(X_train))
  }
  model <- glmnet(X_train, Y_train, alpha = a, family = "cox", lambda = bestlambda, penalty.factor = penalty_factors)
  # Predict on the testing set
  All_test$preds <- predict(model, newx = X_test, s = bestlambda, type = "link")
  
  if (any(is.na(All_test$preds))) {
    cox_formula <- as.formula(Surv(surv, event) ~ 1)
    marker <- 1
    c_index <- NA
  } else {
    cox_formula <- as.formula(Surv(surv, event) ~ preds)
    marker <- All_test$preds
    
    # Calculate c-index using survcomp:
    c_index <- concordance.index(
      x = All_test$preds,
      surv.time = All_test$surv,
      surv.event = All_test$event
    )$c.index
  }
  
  # Calculate AUC using time-dependent ROC
  AUCroc <- timeROC(T = All_test$surv, delta = All_test$event, marker = marker, cause = 1,
                    times = seq(0.9, 4.9, 0.1), ROC = TRUE)$AUC
  # Fit Cox proportional hazards model
  cox_fit <- coxph(cox_formula, data = All_test, x = TRUE, y = TRUE)
  if (any(is.na(coef(cox_fit)))) {
    return(list(
      time_roc = as.numeric(rep(NA, length(seq(0.9, 4.9, 0.1)))),
      ibs = as.numeric(rep(NA, length(seq(0, 4.9, 0.1)))),
      c_index = NA
    ))
  }
  # Calculate Integrated Brier Score using PEC
  pecfit <- pec(object = list(cox_fit), formula = Surv(surv, event) ~ preds, data = All_test,
                exact = FALSE, times = seq(0, 4.9, 0.1), cens.model = "cox", splitMethod = "none", B = 0)
  # Store the Integrated Brier Score (IBS)
  ibs <- crps(pecfit, times = seq(0, 4.9, 0.1), start = 0.1)
  
  return(list(
    time_roc = unname(AUCroc),
    ibs = unname(ibs[1, 1:50]),
    c_index = c_index
  ))
}
