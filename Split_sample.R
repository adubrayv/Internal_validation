##split sample
## source
source("lambda1setuning.R")

## sample division of datasets using 2/3 of training sample and 1/3 of testing sample
###optimization of lambda.1se on training sample and prediction with mean AUC(t), Q10 and Q90 on testing sample

## scenario 1
#lasso-like

AUC_SSV_lasso1 = lapply(1:100, function (i) split_sample_AUC(sim_data100[[i]], 0.95, bestlambda1[[i]], 15000, Weight = FALSE))
result_lasso1 = Get_AUC_SSV(AUC_SSV_lasso1)

#enet

AUC_SSV_enet1 = lapply(1:100, function (i) split_sample_AUC(sim_data100[[i]], 0.5, bestlambdaenet1[[i]], 15000, Weight = FALSE))
result_enet1 = Get_AUC_SSV(AUC_SSV_enet1)

#Ridge-like

AUC_SSV_ridge1 = lapply(1:100, function (i) split_sample_AUC(sim_data100[[i]], 0.05, bestlambdaridge1[[i]], 15000, Weight = FALSE))
result_ridge1 = Get_AUC_SSV(AUC_SSV_ridge1)
result_SSV1 = cbind(result_lasso1, result_enet1, result_ridge1)

##adaptative elastic net

AUC_SSV_adapt1 = lapply(1:100, function (i) split_sample_AUC(sim_data100[[i]], 0.5, bestlambdaenet1[[i]], 15000, Weight = TRUE))
result_SSV_adapt1 = Get_AUC_SSV(AUC_SSV_adapt1)
save(result_SSV_adapt1, file = "result_SSV_adapt1.RData")

## scenario 2
#lasso

AUC_SSV_lasso2 = lapply(1:100, function (i) split_sample_AUC(sim_data500[[i]], 0.95, bestlambda2[[i]], 15000, Weight = FALSE))
result_lasso2 = Get_AUC_SSV(AUC_SSV_lasso2)
save(result_lasso2, file = "result_lasso2.RData")

#enet

AUC_SSV_enet2 = lapply(1:100, function (i) split_sample_AUC(sim_data500[[i]], 0.5, bestlambdaenet2[[i]], 15000, Weight = FALSE))
result_enet2 = Get_AUC_SSV(AUC_SSV_enet2)
save(result_enet2, file = "result_enet2.RData")

#Ridge-like

AUC_SSV_ridge2 = lapply(1:100, function (i) split_sample_AUC(sim_data500[[i]], 0.05, bestlambdaridge2[[i]], 15000, Weight = FALSE))
result_ridge2 = Get_AUC_SSV(AUC_SSV_ridge2)
save(result_ridge2, file = "result_ridge2.RData")
result_SSV2 = cbind(result_lasso2, result_enet2, result_ridge2)

##Adaptative Enet
AUC_SSV_adapt2 = lapply(1:100, function (i) split_sample_AUC(sim_data500[[i]], 0.5, bestlambdaenet2[[i]], 15000, Weight = TRUE))
result_SSV_adapt2 = Get_AUC_SSV(AUC_SSV_adapt2)
save(result_SSV_adapt2, file = "result_SSV_adapt2.RData")

## scenario 3
#lasso-like

AUC_SSV_lasso3 = lapply(1:100, function (i) split_sample_AUC(sim_data1000[[i]], 0.95, bestlambda3[[i]], 15000, Weight = FALSE))
result_lasso3 = Get_AUC_SSV(AUC_SSV_lasso3)
save(result_lasso3, file = "result_lasso3.RData")

#enet

AUC_SSV_enet3 = lapply(1:100, function (i) split_sample_AUC(sim_data1000[[i]], 0.5, bestlambdaenet3[[i]], 15000, Weight = FALSE))
result_enet3 = Get_AUC_SSV(AUC_SSV_enet3)
save(result_enet3, file = "result_enet3.RData")

#Ridge-like

AUC_SSV_ridge3 = lapply(1:100, function (i) split_sample_AUC(sim_data1000[[i]], 0.05, bestlambdaridge3[[i]], 15000, Weight = FALSE))
result_ridge3 = Get_AUC_SSV(AUC_SSV_ridge3)
result_SSV3 = cbind(result_lasso3, result_enet3, result_ridge3)
save(result_SSV3, file = "result_SSV3.RData")

##Adaptative Enet

AUC_SSV_adapt3 = lapply(1:100, function (i) split_sample_AUC(sim_data1000[[i]], 0.5, bestlambdaenet3[[i]], 15000, Weight = TRUE))
result_SSV_adapt3 = Get_AUC_SSV(AUC_SSV_adapt3)
save(result_SSV_adapt3, file = "result_SSV_adapt3.RData")

##saving the results
result_SSV1 = cbind(result_lasso1, result_enet1, result_ridge1, result_SSV_adapt1)
save(result_SSV1, file = "result_SSV1.RData")
result_SSV2 = cbind(result_lasso2, result_enet2, result_ridge2, result_SSV_adapt2)
save(result_SSV2, file = "result_SSV2.RData")
result_SSV3 = cbind(result_lasso3, result_enet3, result_ridge3, result_SSV_adapt3)
save(result_SSV3, file = "result_SSV3.RData")
