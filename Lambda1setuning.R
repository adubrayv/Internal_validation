##model and preparation before bootstrap and cross-validation
##survival object for all scenarios
survobj1 <- lapply(sim_data100, survobj_fun)
survobj2 <- lapply(sim_data500, survobj_fun)
survobj3 <- lapply(sim_data1000, survobj_fun)

##lambda.1se tuning
##lasso-like hypertuning of lambda1se parameters was obtained by cross validation with a = 0.95
##Elastic-Net hypertuning of lambda1se parameters was obtained by cross validation with a = 0.5
##Ridge-like hypertuning of lambda1se parameters was obtained by cross validation with a = 0.05

## lasso-like n = 100
bestlambda1 = lapply(seq_along(sim_data100), function (i) lambda_1se(sim_data100[[i]], 
                                                                     survobj1[[i]], 
                                                                     0.95))
save(bestlambda1, file = "bestlambda1.RData")
fitmodel1 = lapply(seq_along(sim_data100), function (i) fitmodel(sim_data100[[i]],
                                                                 0.95,
                                                                 survobj1[[i]], 
                                                                 bestlambda1[[i]]))

## Enet
bestlambdaenet1 = lapply(seq_along(sim_data100), function (i) lambda_1se(sim_data100[[i]], 
                                                                         survobj1[[i]], 
                                                                         0.5))
save(bestlambdaenet1, file = "bestlambdaenet1.RData")
fitmodelenet1 = lapply(seq_along(sim_data100), function (i) fitmodel(sim_data100[[i]], 
                                                                     0.5,
                                                                     survobj1[[i]], 
                                                                     bestlambdaenet1[[i]]))

## ridge-like
bestlambdaridge1 = lapply(seq_along(sim_data100), function (i) lambda_1se(sim_data100[[i]], 
                                                                          survobj1[[i]], 
                                                                          0.05))
save(bestlambdaridge1, file = "bestlambdaridge1.RData")
fitmodelridge1 = lapply(seq_along(sim_data100), function (i) fitmodel(sim_data100[[i]],
                                                                      0.05,
                                                                      survobj1[[i]], 
                                                                      bestlambdaridge1[[i]]))
##Adaptative Enet
bestlambdaadapt1 = lapply(seq_along(sim_data100), function (i) lambda_1se(sim_data100[[i]], 
                                                                         survobj1[[i]], 
                                                                         0.5))
save(bestlambdaadapt1, file = "bestlambdaadapt1.RData")
fitmodeladapt1 <- lapply(seq_along(sim_data100), function(i) {
  fitmodelw(sim_data100[[i]], 0.5, survobj1[[i]], bestlambdaenet1[[i]])
})
save(fitmodeladapt1, file = "fitmodeladapt1.RData")
## n = 500
##lasso like
bestlambda2 = lapply(seq_along(sim_data500), function (i) lambda_1se(sim_data500[[i]], 
                                                                     survobj2[[i]], 
                                                                     0.95))
save(bestlambda2, file = "bestlambda2.RData")
fitmodel2 = lapply(seq_along(sim_data500), function (i) fitmodel(sim_data500[[i]], 
                                                                 0.95,
                                                                 survobj2[[i]], 
                                                                 bestlambda2[[i]]))

## Enet
bestlambdaenet2 = lapply(seq_along(sim_data500), function (i) lambda_1se(sim_data500[[i]], 
                                                                         survobj2[[i]], 
                                                                         0.5))
save(bestlambdaenet2, file = "bestlambdaenet2.RData")
fitmodelenet2 = lapply(seq_along(sim_data500), function (i) fitmodel(sim_data500[[i]], 
                                                                     0.5,
                                                                     survobj2[[i]], 
                                                                     bestlambdaenet2[[i]]))
##Adaptative Enet
bestlambdaadapt2 = lapply(seq_along(sim_data500), function (i) lambda_1se(sim_data500[[i]], 
                                                                          survobj2[[i]], 
                                                                          0.5))
save(bestlambdaadapt2, file = "bestlambdaadapt2.RData")
fitmodeladapt2 = lapply(seq_along(sim_data500), function (i) fitmodel(sim_data500[[i]],
                                                                      0.5,
                                                                     survobj2[[i]], 
                                                                     bestlambdaadapt2[[i]]))

## ridge-like
bestlambdaridge2 = lapply(seq_along(sim_data500), function (i) lambda_1se(sim_data500[[i]], 
                                                                          survobj2[[i]], 
                                                                          0.05))
save(bestlambdaridge2, file = "bestlambdaridge2.RData")
fitmodelridge2 = lapply(seq_along(sim_data500), function (i) fitmodel(sim_data500[[i]], 
                                                                      0.05,
                                                                      survobj2[[i]], 
                                                                      bestlambdaridge2[[i]]))

## n = 1000
## Lasso-like
bestlambda3 = lapply(seq_along(sim_data1000), function (i) lambda_1se(sim_data1000[[i]], 
                                                                      survobj3[[i]], 
                                                                      0.95))
save(bestlambda3, file = "bestlambda3.RData")
fitmodel3 = lapply(seq_along(sim_data1000), function (i) fitmodel(sim_data1000[[i]],
                                                                  0.95,
                                                                  survobj3[[i]], 
                                                                  bestlambda3[[i]]))

## Enet
bestlambdaenet3 = lapply(seq_along(sim_data1000), function (i) lambda_1se(sim_data1000[[i]], 
                                                                          survobj3[[i]], 
                                                                          0.5))
save(bestlambdaenet3, file = "bestlambdaenet3.RData")
fitmodelenet3 = lapply(seq_along(sim_data1000), function (i) fitmodel(sim_data1000[[i]],
                                                                      0.5,
                                                                      survobj3[[i]], 
                                                                      bestlambdaenet3[[i]]))

## ridge_like
bestlambdaridge3 = lapply(seq_along(sim_data1000), function (i) lambda_1se(sim_data1000[[i]],
                                                                           survobj3[[i]],
                                                                           0.05))
save(bestlambdaridge3, file = "bestlambdaridge3.RData")
fitmodelridge3 = lapply(seq_along(sim_data1000), function (i) fitmodel(sim_data1000[[i]], 
                                                                       0.05,
                                                                       survobj3[[i]], 
                                                                       bestlambdaridge3[[i]]))
##Adaptative Enet
bestlambdaadapt3 = lapply(seq_along(sim_data1000), function (i) lambda_1se(sim_data1000[[i]], 
                                                                          survobj3[[i]], 
                                                                          0.5))
save(bestlambdaadapt3, file = "bestlambdaadapt3.RData")
fitmodeladapt3 = lapply(seq_along(sim_data1000), function (i) fitmodel(sim_data1000[[i]], 
                                                                       0.5,
                                                                     survobj3[[i]], 
                                                                     bestlambdaadapt3[[i]]))

##save the results
save(fitmodel1, file = "fitmodel1.Rdata")
save(fitmodel2, file = "fitmodel2.Rdata")
save(fitmodel3, file = "fitmodel3.Rdata")
save(fitmodelenet1, file = "fitmodelenet1.Rdata")
save(fitmodelenet2, file = "fitmodelenet2.Rdata")
save(fitmodelenet3, file = "fitmodelenet3.Rdata")
save(fitmodelridge1, file = "fitmodelridge1.Rdata")
save(fitmodelridge2, file = "fitmodelridge2.Rdata")
save(fitmodelridge3, file = "fitmodelridge3.Rdata")