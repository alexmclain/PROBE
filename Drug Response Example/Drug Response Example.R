# unhidem application with CCLE data for all methods except conformal
# Jackknife, which was ran separately due it's computational demand.
#
# This function runs 6 analysis methods (UNHIDEM, MCP, SCAD, LASSO
# ALASSO (adaptive LASSO) and LASSO with the conformal split method),
# 8 drugs and 10 folds (all together 480 analyses).  It takes about
# 2.5-3 hours to run. The results have been saved in the file:
# "Drug_results.rds".
# The 10 folds are used to estimate the out-of-sample prediction
# error of the methods, and the out-of-sample empirical coverage 
# probabilities of 95% prediction intervals.

rm(list=ls())
require(caret)
library(unhidem)
library(ncvreg)
library(sparsevb)
library(varbvs)
library(SSLASSO)
library(plyr)
library(dplyr)
library(glmnet)
library(writexl)
library(conformalInference)

#Reading in data
z_data <- read.csv("CCLE_Zdata.csv")

y_data <- read.csv("CCLE_Ydata.csv")
y_data <- y_data[, colSums(apply(y_data, 2, is.na))==0]

#Prepare objects to store results
nfolds <- 10
mspe = mad <- array(numeric(), dim = c(nfolds, ncol(y_data), 9))
pred_data <- array(numeric(), dim = c(nrow(y_data), ncol(y_data), 9)) 
test_data <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 
ecp       <- array(numeric(), dim = c(nrow(y_data), ncol(y_data), 2))
PI_len    <- array(numeric(), dim = c(nrow(y_data), ncol(y_data), 2)) 

#Run Models
alpha <- 0.05
for(i in 1:ncol(y_data)) {
  
  outcome <- y_data[,i]
  cat(paste0("Drug ", i," of ",ncol(y_data),": ", names(y_data)[i]), "\n")
  
  t_pred_data = t_test_data = t_ecp = t_len = t_ecp_lass = t_ecp_rid = 
    t_len_lass = t_len_rid <- NULL
  
  #Create 10 folds for cross-validation
  set.seed(123)
  ind <- sample(1:(length(outcome)) %% nfolds)
  ind[ind == 0] <- nfolds
  for(j in 1:nfolds){
    
    cat(paste0("Fitting fold ", j," for: "))
    
    #UNHIDEM SECTION
    cat("UNHIDEM, ")
    Z <- as.matrix(z_data[ind!=j,])
    M <- ncol(Z)
    Z_mean <- mean(array(Z))
    Y <- outcome[ind!=j]
    
    Z_test <- as.matrix(z_data[ind==j,])
    Y_test <- outcome[ind==j]
    t_test_data <- c(t_test_data, Y_test)
    
    Z <- Z - Z_mean
    Z_test <- Z_test - Z_mean
    
    results <- unhidem(Y = Y, Z = Z, alpha = alpha) 
    
    pred_res_test <- predict_unhidem_func(results, Z = Z_test, alpha = alpha)
    yhat_unhidem <- pred_res_test$Pred
    
    mspe[j,i,1] <- mean( (Y_test - yhat_unhidem)^2 )
    mad[j,i,1] <- median( abs(Y_test - yhat_unhidem) )
    
    # Empirical coverage probabilities and PI lengths
    t_ecp <- c( t_ecp, 1*I(Y_test > pred_res_test$PI_L & 
                             Y_test < pred_res_test$PI_U) )
    t_len <- c( t_len, pred_res_test$PI_U - pred_res_test$PI_L)
    
    if(FALSE){
    #MCP AND SCAD SECTION
    cat("MCP-SCAD, ")
    X <- Z
    X_test <- Z_test
    
    set.seed(123)
    cv.out_mcp <- cv.ncvreg(X, Y, nfolds = 10,
                            family = "gaussian",penalty = "MCP",seed=123)
    set.seed(123)
    cv.out_scad <- cv.ncvreg(X, Y, nfolds = 10,
                             family = "gaussian",penalty = "SCAD",seed=123)
    
    yhat_mcp  <- predict(cv.out_mcp, X_test, which = cv.out_mcp$min)
    yhat_scad <- predict(cv.out_scad, X_test, which = cv.out_scad$min)
    
    mspe[j,i,2] <- mean( (Y_test - yhat_mcp)^2 )
    mspe[j,i,3] <- mean( (Y_test - yhat_scad)^2 )
    
    mad[j,i,2] <- median( abs(Y_test - yhat_mcp) )
    mad[j,i,3] <- median( abs(Y_test - yhat_scad) )
    
    #LASSO AND ALASSO SECTION
    cat("LASSO-ALASSO, ")
    set.seed(123)
    ridge1_cv <- cv.glmnet(X, Y, alpha = 0, nfolds = 10)
    best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
    set.seed(123)
    cv.out_alasso=cv.glmnet(X, Y, alpha = 1, nfolds = 10,
                            penalty.factor = 1 / abs(best_ridge_coef), keep = TRUE)
    
    set.seed(123)    
    cv.out_lasso=cv.glmnet(X, Y, alpha = 1, nfolds = 10)
    
    yhat_alasso  <- predict(cv.out_alasso, newx = X_test, s = "lambda.min")
    yhat_lasso   <- predict( cv.out_lasso, newx = X_test, s = "lambda.min")
    
    mspe[j,i,4] <- mean( (Y_test - yhat_lasso)^2 )
    mspe[j,i,5] <- mean( (Y_test - yhat_alasso)^2 )
    
    mad[j,i,4] <- median( abs(Y_test - yhat_lasso) )
    mad[j,i,5] <- median( abs(Y_test - yhat_alasso) )
    
    lasso_coefs <- coef(cv.out_lasso,s="lambda.min")[-1]
    lasso_M1 <- length(lasso_coefs[lasso_coefs!=0])
    
    #LASSO SPLIT CONFORMAL SECTION
    cat("Conformal Split, ")
    
    funs = lasso.funs(nlambda=100,cv=TRUE, lambda.min.ratio = 1e-04)
    
    out.split.las <- conformal.pred.split(Z, Y, Z_test, alpha = alpha, seed = 123,
                                          train.fun = funs$train, predict.fun = 
                                            funs$predict, verb = FALSE)
    
    t_ecp_lass <- c( t_ecp_lass, 1*I(Y_test > out.split.las$lo & Y_test < 
                                       out.split.las$up) )
    t_len_lass <- c( t_len_lass, out.split.las$up - out.split.las$lo)
    
    yhat_lasso_conf  <- out.split.las$pred
    
    mspe[j,i,6] <- mean( (Y_test - yhat_lasso_conf)^2 )
    mad[j,i,6] <- median( abs(Y_test - yhat_lasso_conf) )
    
    ### VARBVS
    cat("varbvs, ")
    mod.out <- varbvs(X=Z, Z=NULL, y=Y, family = "gaussian", verbose = FALSE)
    yhat_varbvs <- predict(mod.out,X = Z_test)
    
    mspe[j,i,7] <- mean( (Y_test - yhat_varbvs)^2 )
    mad[j,i,7] <- median( abs(Y_test - yhat_varbvs) )
    
    #### SSLASSO
    cat("SSLASSO, ")
    L <- 400
    mod.out <- SSLASSO(X=Z, y=Y, variance = "unknown", lambda1 = 0.01, 
                                            lambda0 = seq(0.01,M,length.out=L), 
                                            a = lasso_M1, b = M-lasso_M1)
    SSLASSO_coefs <- mod.out$beta[,L]
    yhat_sslasso<- Z_test%*%c(SSLASSO_coefs) + mod.out$intercept[L]
    
    mspe[j,i,8] <- mean( (Y_test - yhat_sslasso)^2 )
    mad[j,i,8] <- median( abs(Y_test - yhat_sslasso) )
    
    ### sparsevb
    cat("and sparsevb. \n")
    mod.out <- svb.fit(X=Z, Y=Y, family = "linear", slab = "laplace", mu = lasso_coefs,
                       intercept = TRUE, alpha = lasso_M1, beta = M - lasso_M1)
    sparsevb_coefs <- mod.out$mu * mod.out$gamma #approximate posterior mean
    yhat_sparsevb<- Z_test%*%c(sparsevb_coefs) + mod.out$intercept
    
    mspe[j,i,9] <- mean( (Y_test - yhat_sparsevb)^2 )
    mad[j,i,9] <- median( abs(Y_test - yhat_sparsevb) )
    }
    
    yhat_mcp <- yhat_unhidem
    yhat_scad <- yhat_unhidem
    yhat_lasso <- yhat_unhidem
    yhat_alasso <- yhat_unhidem
    yhat_lasso_conf <- yhat_unhidem
    yhat_varbvs <- yhat_unhidem
    yhat_sslasso <- yhat_unhidem
    yhat_sparsevb <- yhat_unhidem
     
    #Gather all of the predicted values
    t_pred_data <- rbind( t_pred_data, cbind(yhat_unhidem, yhat_mcp, 
                                             yhat_scad, yhat_lasso, yhat_alasso, 
                                             yhat_lasso_conf, yhat_varbvs, yhat_sslasso, 
                                             yhat_sparsevb) )
  }
  
  pred_data[,i,] <- t_pred_data 
  test_data[,i]  <-  t_test_data
  ecp[,i,1]      <-  t_ecp
  ecp[,i,2]      <-  t_ecp_lass 
  PI_len[,i,1]   <-  t_len 
  PI_len[,i,2]   <-  t_len_lass 
  
}


results <- list("pred_data" = data.frame(pred_data), "test_data" = 
                  data.frame(test_data), "ecp" = data.frame(ecp), 
                "PI_len" =data.frame(PI_len), "mspe" = data.frame(mspe), 
                "mad" = data.frame(mad))

saveRDS(results,"Drug_results_W2.rds")


