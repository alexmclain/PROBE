
setwd("C:/Users/mclaina/OneDrive - University of South Carolina/Research/Imaging/PROBE/Programs/PROBE/Drug Response Example")

library("conformalInference")
require(caret)
library(ncvreg)
library(dplyr)
library(glmnet)
library(writexl)

#Reading in data
z_data <- read.csv("CCLE_Xdata.csv")

y_data <- read.csv("CCLE_Ydata.csv")
y_data <- y_data[, colSums(apply(y_data, 2, is.na))==0]


nfolds <- 10
mse = mad <- array(numeric(), dim = c(nfolds, ncol(y_data)))
pred_data <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 
ecp = PI_len  <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 

test_data <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 

#Run Models
for(i in 1:ncol(y_data)) {
  
  outcome <- y_data[,i]
  cat(paste0("Drug: ", names(y_data)[i]), "\n")
  
  t_pred_data = t_test_data = t_ecp_lass = t_len_lass  <- NULL
  
  #Create 10 folds for cross-validation
  set.seed(123)
  ind <- sample(1:(length(outcome)) %% nfolds)
  ind[ind == 0] <- nfolds
  for(j in 1:nfolds){
    
    print(paste0("Fold: ", j))
    
    Z <- as.matrix(z_data[ind!=j,])
    Y <- outcome[ind!=j]
    
    Z_test <- as.matrix(z_data[ind==j,])
    Y_test <- outcome[ind==j]
    t_test_data <- c(t_test_data, Y_test)
    
    
    #LASSO AND RIDGE CONFORMAL SECTION
    print("Lasso Conformal Section")
    
    
    funs = lasso.funs(cv = TRUE)
    set.seed(123)
    system.time(out.jack.las <- conformal.pred.jack(Z, Y, Z_test, alpha=0.05,
                                          train.fun=funs$train, predict.fun=funs$predict, 
                                          verb = TRUE))
    
    t_ecp_lass <- c( t_ecp_lass, 1*I(Y_test > out.jack.las$lo & Y_test < out.jack.las$up ) )
    t_len_lass <- c( t_len_lass, out.jack.las$up - out.jack.las$lo)
    
    yhat_lasso_conf  <- out.jack.las$pred
    
    mse[j,i] <- mean( (Y_test - yhat_lasso_conf)^2 )
    mad[j,i] <- median( abs(Y_test - yhat_lasso_conf) )
    
    #Gather all of the predicted values
    t_pred_data <- rbind( t_pred_data, yhat_lasso_conf)
  }
  
  pred_data[,i] <-  t_pred_data 
  test_data[,i] <-  t_test_data
  ecp[,i]       <-  t_ecp_lass 
  PI_len[,i]    <-  t_len_lass 
  
}


results <- list("Pred" = data.frame(pred_data), "test" = data.frame(test_data), 
               "ecp" = data.frame(ecp), "PI_Length" =data.frame(PI_len))

saveRDS(results,"Conform_lasso_jack.rds")





