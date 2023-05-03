
setwd("C:/Users/mclaina/OneDrive - University of South Carolina/Research/Imaging/PROBE/Programs/PROBE/Drug Response Example")


library(ebreg)
require(caret)
library(dplyr)
library(writexl)

#Reading in data
z_data <- read.csv("CCLE_Xdata.csv")

y_data <- read.csv("CCLE_Ydata.csv")
y_data <- y_data[, colSums(apply(y_data, 2, is.na))==0]


nfolds <- 10
p <- ncol(z_data)
mspe = mad <- array(numeric(), dim = c(nfolds, ncol(y_data)))
pred_data <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 
ecp = PI_len  <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 

test_data <- array(numeric(), dim = c(nrow(y_data), ncol(y_data))) 

#Run Models
for(i in 1:ncol(y_data)) {
  
  outcome <- y_data[,i]
  cat(paste0("Drug: ", names(y_data)[i]), "\n")
  
  t_pred_data = t_test_data = t_ecp_ebreg = t_len_ebreg <- NULL
  
  #Create 10 folds for cross-validation
  set.seed(123)
  ind <- sample(1:(length(outcome)) %% nfolds)
  ind[ind == 0] <- nfolds
  for(j in 1:nfolds){
    
    cat("Fold",j,"\n")
    
    Z <- as.matrix(z_data[ind!=j,])
    Y <- outcome[ind!=j]
    
    z_mean <- apply(Z, 2, mean)
    y_mean <- mean(Y)
    
    Z_test <- as.matrix(z_data[ind==j,])
    Y_test <- outcome[ind==j]
    t_test_data <- c(t_test_data, Y_test)
    
    
    #empirical Bayes method
    n <- length(Y)
    log.f <- function(x) log(1/n) + log(x <= n) #log prior of the model size
    
    system.time(out.ebreg <- ebreg( Y, Z, Z_test, standardized = FALSE, alpha = 0.95, sig2=NULL, prior = TRUE, 
                   log.f = log.f, M=5000, sample.beta = FALSE, pred = TRUE, 
                   conf.level = 0.95))
    
    
    yhat_ebreg  <- out.ebreg$ynew.mean
    
    mspe[j,i] <- mean( (Y_test - yhat_ebreg)^2 )
    mad[j,i] <- median( abs(Y_test - yhat_ebreg) )
    
    t_pred_data <- c( t_pred_data, yhat_ebreg)
    t_ecp_ebreg <- c( t_ecp_ebreg, 1*I(Y_test > out.ebreg$PI[1,] & Y_test < out.ebreg$PI[2,]) )
    t_len_ebreg <- c( t_len_ebreg, out.ebreg$PI[2,] - out.ebreg$PI[1,])
    
    t_res <- c(mean( (Y_test - yhat_ebreg)^2 ), median( abs(Y_test - yhat_ebreg) ), 
               mean(1*I(Y_test > out.ebreg$PI[1,] & Y_test < out.ebreg$PI[2,])))
    cat(t_res,"\n")
  }
  
  pred_data[,i] <-  t_pred_data 
  test_data[,i] <-  t_test_data
  ecp[,i]       <-  t_ecp_ebreg 
  PI_len[,i]    <-  t_len_ebreg 
  
}


results <- list("Pred" = data.frame(pred_data), "test" = data.frame(test_data), 
               "ecp" = data.frame(ecp), "PI_Length" =data.frame(PI_len), "mspe" = mspe, "mad" = mad)

saveRDS(results,"EBreg_unif.rds")





