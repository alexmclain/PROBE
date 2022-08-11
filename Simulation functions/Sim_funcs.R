
library(RandomFields)



simulation_func <- function(args, parlist, B, bin, ebreg_I, verbose = TRUE){
  args_list <- as.numeric(unlist(parlist[args,]))
  
  NUM <- as.numeric(args_list[1])
  x <- seq(1, NUM, 1)
  K <- B
  N <- 400
  M <- NUM^2
  M1 <- as.integer(M*as.numeric(args_list[2]))
  p <- 1-as.numeric(args_list[2])
  err_sd <- 0.75
  alpha <- 0.05
  eta <-   as.numeric(args_list[3]) #Eta value
  eta_var <- (eta*sqrt(pi)/sqrt(2))^2
  sig_nois <- as.numeric(args_list[4])
  lat_sp <- 10
  sig_sp <- 20
  lat_data <- 0
  RFoptions(spConform=FALSE)
  
  ### Generating a big dataset to get the variance of the signal. First run is required so
  ### the seed works appropriately (bug in RandomFields package).
  data <- sim_bin_LR_2(lat_data,err_sd, x, 1,M,M1,sig_sp,lat_sp,eta_var,seed = 238476 )
  data <- sim_bin_LR_2(lat_data,err_sd, x, 10000,M,M1,sig_sp,lat_sp,eta_var,seed = 238476 )
  #Signal and beta coefficients
  set.seed(args)
  signal <- data$signal
  sig_ind <- data$sig_ind
  t_eta <-  runif(M, 0, 2*eta)# data$eta#
  
  #Z data
  if(bin){Z <- matrix(t(data$Z), nrow = N, ncol = M, byrow = TRUE)}else{Z <- data$Z_cont}
  
  # Generate outcome
  eta_vec <- sig_ind*array(t_eta)
  eta_i <- apply(t(Z)*c(eta_vec),2,sum)
  sig_var <- var(eta_i)
  
  sigma <- (sig_var/sig_nois)^(0.5)
  
  
  
  UNHIDEM_res <- matrix(0,K,13)
  colnames(UNHIDEM_res) <- c("Sum_delta", "Sum_delta_sig", "MSE", "MAD", "ECP_PI", 
                             "ECP_CI", "Iter", "Sigma2_est", "Test_MSPE", "Test_MSE", "Beta_err", 
                             "Conv", "Time")
  
  varbvs_res <- SSLASSO_res <- sparsevb_res <- sparsevb_c_res <- matrix(-999,K,8)
  LASSO_res <- adap_LASSO_res <- SCAD_res <- MCP_res <- matrix(-999,K,8)
  ebreg_res <- matrix(-999,K,9)
  
  colnames(varbvs_res) <- c("Varbvs_Signals", "Varbvs_Corr_Signals", "Varbvs_MSE", "Varbvs_MAD", 
                            "Varbvs_Test_MSPE", "Varbvs_Test_MSE", "Varbvs_Beta_err", "Varbvs_time")
  
  colnames(sparsevb_res) <- c("sparsevb_Signals", "sparsevb_Corr_Signals", "sparsevb_MSE", 
                              "sparsevb_MAD", "sparsevb_Test_MSPE", "sparsevb_Test_MSE", 
                              "sparsevb_Beta_err", "sparsevb_time")
  
  colnames(sparsevb_c_res) <- c("sparsevb_c_Signals", "sparsevb_c_Corr_Signals", "sparsevb_c_MSE", 
                                "sparsevb_c_MAD", "sparsevb_c_Test_MSPE", "sparsevb_c_Test_MSE", 
                                "sparsevb_c_Beta_err", "sparsevb_c_time")
  
  colnames(SSLASSO_res) <- c("SSLASSO_Signals", "SSLASSO_Corr_Signals", "SSLASSO_MSE", "SSLASSO_MAD", 
                             "SSLASSO_Test_MSPE", "SSLASSO_Test_MSE", "SSLASSO_Beta_err", 
                             "SSLASSO_time")
  
  colnames(ebreg_res) <- c("ebreg_Signals", "ebreg_Corr_Signals", "ebreg_res_MSE", "ebreg_MAD", 
                           "ebreg_Test_MSPE", "ebreg_Test_MSE", "ebreg_Beta_err", "ebreg_ECP_PI",
                           "ebreg_time")
  
  colnames(LASSO_res) <- c("lasso_Sum_delta","lasso_Sum_delta_sig", 
                           "lasso_MSE", "lasso_MAD", 
                           "lasso_Obs_test_MSPE", "lasso_test_MSE", "lasso_Beta_err", "lasso_time")
  
  colnames(adap_LASSO_res) <- c("adap_lasso_Sum_delta","adap_lasso_Sum_delta_sig", 
                                "adap_lasso_MSE", "adap_lasso_MAD", 
                                "adap_lasso_Obs_test_MSPE", "adap_lasso_test_MSE", 
                                "adap_lasso_Beta_err", "adap_lasso_time")
  
  colnames(SCAD_res) <- c("SCAD_Sum_delta","SCAD_Sum_delta_sig", 
                          "SCAD_MSE", "SCAD_MAD", 
                          "SCAD_Obs_test_MSPE", "SCAD_test_MSE", "SCAD_Beta_err", "SCAD_time")
  colnames(MCP_res) <- c("MCP_Sum_delta","MCP_Sum_delta_sig", 
                         "MCP_MSE", "MCP_MAD", 
                         "MCP_Obs_test_MSPE", "MCP_test_MSE", "MCP_Beta_err", "MCP_time")
  
  if(!ebreg_I){ebreg_res <- NULL}
  start <- 1000*(args-1)
  ## Set convergence criteria
  maxit <- 1500
  ep <- 0.1
  
  
  if(bin){bin_text <- c("bin")}else{bin_text <- c("cont")}
  filname <- paste0("All_sim",bin_text," eta ",eta," M ",M," SNR ", 
                    sig_nois," M1 ",M1,".csv")
  print(filname)
  
  for(k in 1:K){
    
    # Generate Z and Signal Data
    set.seed(start + k)
    data <- sim_bin_LR_2(lat_data,err_sd,x,N,M,M1,sig_sp,lat_sp,eta_var,seed = start + k)
    
    #Signal and beta coefficients
    signal <- data$signal
    sig_ind <- data$sig_ind
    t_eta <-  runif(M, 0, 2*eta)# data$eta#
    LP_data <- expand.grid(x1=x,x2=x)
    LP_data$signal <- sig_ind
    LP_data$Signal <- factor((LP_data$signal-1))
    
    #Z data
    if(bin){Z <- matrix(t(data$Z), nrow = N, ncol = M, byrow = TRUE)}else{Z <- data$Z_cont}
    
    # Generate outcome
    eta_vec <- sig_ind*array(t_eta)
    eta_i <- apply(t(Z)*c(eta_vec),2,sum)
    Y <- eta_i + rnorm(N,0,sigma)
    
    # Generate Test Data
    set.seed(start + k + 12321)
    t_data <- RFsimulate(model = RMstable(alpha = 2, scale = lat_sp,var = 1), x=x, y=x,grid=TRUE,n=N)
    t_datamat <- t(matrix(array(t_data),M,N))+rnorm(N)*err_sd
    if(bin){
      Z_test <- matrix(1*I(t(t(t_datamat)+c(array(lat_data))) < 0), N, M)
    }else{
      Z_test <- t(t(t_datamat)+c(array(lat_data))) 
    }
    eta_test <- apply(t(Z_test)*c(eta_vec),2,sum)
    Y_test <- eta_test + rnorm(N,0,sigma)
    
    test3 <- system.time(mod.out <- unhidem(Y = Y, Z = Z, alpha = alpha, eta_i = eta_i, 
                                            signal = signal))
    
    alpha_est <- mod.out$Calb_mod$coef[2]
    gamma_est <- mod.out$E_step$delta
    beta_est <- mod.out$beta_hat
    beta_ast_est <- gamma_est*alpha_est*beta_est
    sigma2_est_new <- mod.out$Calb_mod$sigma2_est
    Y_pred <- mod.out$Calb_mod$Y_pred
    
    #Prediction for test data
    pred_res_test <- predict_unhidem_func(mod.out, Z_test, X = NULL, alpha = alpha)
    # Proportion of test PIs that contain the test observation
    ECP_PI <- mean(1*I(Y_test>pred_res_test$PI_L & Y_test<pred_res_test$PI_U))
    # Proportion of test CIs that contain the test true signal
    ECP_CI <- mean(1*I(eta_test>pred_res_test$CI_L & eta_test<pred_res_test$CI_U))
    
    UNHIDEM_res[k,] <- c(sum(gamma_est),
                         sum(gamma_est[signal]), 
                         mean((Y_pred-eta_i)^2), 
                         median(abs(Y_pred - eta_i)), ECP_PI, ECP_CI,   
                         mod.out$count, sigma2_est_new, mean((Y_test - pred_res_test$Pred)^2), 
                         mean((eta_test - pred_res_test$Pred)^2), 
                         mean((eta_vec - beta_ast_est)^2),
                         mod.out$conv, test3[3])
    
    
    #### LASSO
    test3 <- system.time(cv.out <- cv.glmnet(Z,Y,alpha = 1,nfolds = 10,lambda.min.ratio=0.001))
    while(cv.out$lambda.min==min(cv.out$lambda)){
      ### Further decreasing the lambda if the smallest lambda was the best.
      cv.out <- cv.glmnet(Z,Y,alpha = 1,nfolds = 10,lambda=exp(c(log(cv.out$lambda.1se),seq(log(min(cv.out$lambda)),log(min(cv.out$lambda)*0.001),length.out = 100))))
    }
    
    lasso_coefs <- coef(cv.out,s="lambda.min")[-1]
    lasso_pred  <- predict(cv.out,newx = Z,s="lambda.min")
    lasso_mse <- mean((lasso_pred - eta_i)^2)
    lasso_pred_test<- predict(cv.out,newx = Z_test,s="lambda.min")
    lasso_mad <- median(abs(lasso_pred_test - eta_test))
    lasso_mspe     <- mean((lasso_pred_test - Y_test)^2)
    lasso_mse_test <- mean((lasso_pred_test - eta_test)^2)
    
    
    LASSO_res[k,] <- c(length(lasso_coefs[lasso_coefs!=0]),
                       length(lasso_coefs[lasso_coefs!=0 & sig_ind==1]), 
                       lasso_mse, lasso_mad , lasso_mspe, lasso_mse_test, 
                       mean((lasso_coefs - eta_vec)^2), test3[3])
    
    #### ADAPTIVE LASSO
    test3 <- system.time(t1 <- try(cv.out <- adap_lasso(Y,Z), silent = TRUE))
    
    if(!(attr(t1,"class")=="try-error")){
      lasso_coefs <- coef(cv.out,s="lambda.min")[-1]
      lasso_pred  <- predict(cv.out,newx = Z,s="lambda.min")
      lasso_mse <- mean((lasso_pred - eta_i)^2)
      lasso_mad <- median(abs(lasso_pred - eta_i))
      lasso_pred_test<- predict(cv.out,newx = Z_test,s="lambda.min")
      lasso_mspe     <- mean((lasso_pred_test - Y_test)^2)
      lasso_mse_test <- mean((lasso_pred_test - eta_test)^2)
      
      
      adap_LASSO_res[k,] <- c(length(lasso_coefs[lasso_coefs!=0]),
                              length(lasso_coefs[lasso_coefs!=0 & sig_ind==1]), 
                              lasso_mse, lasso_mad , lasso_mspe, lasso_mse_test, 
                              mean((lasso_coefs - eta_vec)^2), test3[3])
    }
    
    test3 <- system.time(t1 <- try(cv.out2 <- scad_func(Y,Z,N,M), silent = TRUE))
    
    if(!(attr(t1,"class")=="try-error")){
      scad_coefs <- cv.out2$fit$beta[,cv.out2$min][-1]
      ash_scad <- predict(cv.out2,matrix(as.numeric(Z),N,M),which = cv.out2$min)
      scad_mse <- mean((ash_scad - eta_i)^2)
      scad_mad <- median(abs(ash_scad - eta_i))
      ash_scad_test <- predict(cv.out2,matrix(as.numeric(Z_test),N,M),which = cv.out2$min)
      scad_mspe     <- mean((ash_scad_test - Y_test)^2)
      scad_mse_test <- mean((ash_scad_test - eta_test)^2)
      
      SCAD_res[k,] <- c(length(scad_coefs[scad_coefs!=0]),
                        length(scad_coefs[scad_coefs!=0 & sig_ind==1]), 
                        scad_mse, scad_mad , scad_mspe, scad_mse_test, 
                        mean((scad_coefs - eta_vec)^2), test3[3])
    }
    
    test3 <- system.time(t1 <- try(cv.out3 <- MCP_func(Y,Z,N,M), silent = TRUE))
    
    if(!(attr(t1,"class")=="try-error")){
      mcp_coefs <- cv.out3$fit$beta[,cv.out3$min][-1]
      ash_mcp <- predict(cv.out3,matrix(as.numeric(Z),N,M),which = cv.out3$min)
      mcp_mse <- mean((ash_mcp - eta_i)^2)
      ash_mcp_test <- predict(cv.out3,matrix(as.numeric(Z_test),N,M),which = cv.out3$min)
      mcp_mspe     <- mean((ash_mcp_test - Y_test)^2)
      mcp_mse_test <- mean((ash_mcp_test - eta_test)^2)
      mcp_mad      <- median(abs(ash_mcp_test - eta_test))
      
      MCP_res[k,] <- c(length(mcp_coefs[mcp_coefs!=0]),
                       length(mcp_coefs[mcp_coefs!=0 & sig_ind==1]), 
                       mcp_mse, mcp_mad , mcp_mspe, mcp_mse_test, 
                       mean((mcp_coefs - eta_vec)^2), test3[3])
    }
    
    #### varbvs
    test3 <- system.time(mod.out <- varbvs(X=Z, Z=NULL, y=Y, 
                                           family = "gaussian", verbose = FALSE))
    ## 
    varbvs_pip   <- mod.out$pip
    varbvs_coefs <- mod.out$beta*varbvs_pip
    varbvs_pred  <- predict(mod.out,X = Z)
    varbvs_mse <- mean((varbvs_pred - eta_i)^2)
    varbvs_mad <- median(abs(varbvs_pred - eta_i))
    varbvs_pred_test<- predict(mod.out,X = Z_test)
    varbvs_mspe     <- mean((varbvs_pred_test - Y_test)^2)
    varbvs_mse_test <- mean((varbvs_pred_test - eta_test)^2)
    
    
    varbvs_res[k,] <- c(sum(varbvs_pip), sum(varbvs_pip*sig_ind), 
                        varbvs_mse, varbvs_mad , varbvs_mspe, varbvs_mse_test, 
                        mean((varbvs_coefs - eta_vec)^2), test3[3])
    
    #### SSLASSO
    L <- 400
    test3 <- system.time(t1 <- try(mod.out <- SSLASSO(X=Z, y=Y, variance = "unknown", lambda1 = 0.01, 
                                                      lambda0 = seq(0.01,M,length.out=L), 
                                                      a = M1, b = M-M1), silent = TRUE))
    if(attr(t1,"class") == "SSLASSO"){
      SSLASSO_coefs <- mod.out$beta[,L]
      SSLASSO_pred  <- Z%*%c(SSLASSO_coefs) + mod.out$intercept[L]
      SSLASSO_mse <- mean((SSLASSO_pred - eta_i)^2)
      SSLASSO_mad <- median(abs(SSLASSO_pred - eta_i))
      SSLASSO_pred_test<- Z_test%*%c(SSLASSO_coefs) + mod.out$intercept[L]
      SSLASSO_mspe     <- mean((SSLASSO_pred_test - Y_test)^2)
      SSLASSO_mse_test <- mean((SSLASSO_pred_test - eta_test)^2)
      
      
      SSLASSO_res[k,] <- c(length(mod.out$model),
                           length(mod.out$model %in% signal), 
                           SSLASSO_mse, SSLASSO_mad , SSLASSO_mspe, SSLASSO_mse_test, 
                           mean((SSLASSO_coefs - eta_vec)^2), test3[3])
    }
    
    #empirical Bayes method
    if(ebreg_I){
      n <- length(Y)
      log.f <- function(x) log(1/n) + log(x <= n) #log prior of the model size
      
      test3 <- system.time(t1 <- try(out.ebreg <- ebreg( Y, Z, rbind(Z,Z_test), standardized = FALSE, 
                                                         alpha = 0.99, prior = TRUE, M = 2000, 
                                                         log.f = log.f, sample.beta = TRUE, 
                                                         pred = TRUE, conf.level = 0.95), silent = TRUE))
      
      if(is.null(attr(t1,"class"))){
        
        yhat_ebreg  <- out.ebreg$ynew.mean
        # Proportion of test PIs that contain the test observation
        ECP_PI <- mean(1*I(Y_test> out.ebreg$PI[1,401:800] & Y_test < out.ebreg$PI[2,401:800]))
        
        ebreg_pip   <- out.ebreg$incl.prob
        ebreg_coefs <- out.ebreg$beta.mean*ebreg_pip
        ebreg_pred  <- yhat_ebreg[1:400]
        
        ebreg_mse <- mean((ebreg_pred - eta_i)^2)
        ebreg_mad <- median(abs(ebreg_pred - eta_i))
        
        ebreg_pred_test  <- yhat_ebreg[401:800]
        ebreg_mspe     <- mean((ebreg_pred_test - Y_test)^2 )
        ebreg_mse_test <- mean((ebreg_pred_test - eta_test)^2)
        
        ebreg_res[k,] <- c(sum(ebreg_pip),
                           sum(ebreg_pip*sig_ind), 
                           ebreg_mse, ebreg_mad , ebreg_mspe, ebreg_mse_test, 
                           mean((ebreg_coefs - eta_vec)^2), ECP_PI, test3[3])
      }
    }
    
    
    test3 <- system.time(t1 <- try(test <- svb.fit(X=Z, Y=Y, family = "linear", slab = "laplace",
                                                   intercept = TRUE), silent = TRUE))
    
    if(is.null(attr(t1,"class"))){
      sparsevb_coefs <- test$mu * test$gamma #approximate posterior mean
      sparsevb_pred  <- Z%*%c(sparsevb_coefs) + test$intercept
      sparsevb_mse <- mean((sparsevb_pred - eta_i)^2)
      sparsevb_mad <- median(abs(sparsevb_pred - eta_i))
      sparsevb_pred_test<- Z_test%*%c(sparsevb_coefs) + test$intercept
      sparsevb_mspe     <- mean((sparsevb_pred_test - Y_test)^2)
      sparsevb_mse_test <- mean((sparsevb_pred_test - eta_test)^2)
      
      
      sparsevb_res[k,] <- c(sum(test$gamma ),
                            sum(test$gamma*sig_ind), 
                            sparsevb_mse, sparsevb_mad , sparsevb_mspe, sparsevb_mse_test, 
                            mean((sparsevb_coefs - eta_vec)^2), test3[3])
      
      sparsevb_c_res[k,] <- sparsevb_res[k,]
      
      if(max(sparsevb_coefs)>0){
        sparsevb_pred_cov  <- Z%*%c(sparsevb_coefs) + test$intercept
        sparse_mod <- lm(Y~sparsevb_pred_cov)
        sparsevb_coefs <- sparsevb_coefs*sparse_mod$coefficients[2] #approximate posterior mean
        sparsevb_pred <- sparse_mod$fitted.values
        
        sparsevb_mse <- mean((sparsevb_pred - eta_i)^2)
        sparsevb_mad <- median(abs(sparsevb_pred - eta_i))
        sparsevb_pred_test <- Z_test%*%c(sparsevb_coefs) + 
          sparse_mod$coefficients[2]*test$intercept + 
          sparse_mod$coefficients[1]
        sparsevb_mspe     <- mean((sparsevb_pred_test - Y_test)^2)
        sparsevb_mse_test <- mean((sparsevb_pred_test - eta_test)^2)
        
        
        sparsevb_c_res[k,] <- c(sum(test$gamma ),
                                sum(test$gamma*sig_ind), 
                                sparsevb_mse, sparsevb_mad , sparsevb_mspe, sparsevb_mse_test, 
                                mean((sparsevb_coefs - eta_vec)^2), test3[3])
      }
    }
    
    if(ebreg_I){
      k_res <- data.frame(rbind(round(UNHIDEM_res[k,-c(5:8,12)],3),round(LASSO_res[k,],3),
                                round(adap_LASSO_res[k,],3),round(SCAD_res[k,],3),
                                round(MCP_res[k,],3),round(SSLASSO_res[k,],3),
                                round(varbvs_res[k,],3),round(sparsevb_res[k,],3),
                                round(sparsevb_c_res[k,],3),round(ebreg_res[k,-8],3)))
      colnames(k_res)[1:2] <- c("Total_Beta", "Correct_Beta")
      rownames(k_res) <- c("UNHIDEM","LASSO","ALASSO","SCAD","MCP","SSLASSO", 
                           "VARBVS","SPARSEVB", "SPARSEVB_c","EBREG")}
    if(!ebreg_I){
      k_res <- data.frame(rbind(round(UNHIDEM_res[k,-c(5:8,12)],3),round(LASSO_res[k,],3),
                                round(adap_LASSO_res[k,],3),round(SCAD_res[k,],3),
                                round(MCP_res[k,],3),round(SSLASSO_res[k,],3),
                                round(varbvs_res[k,],3),round(sparsevb_res[k,],3),
                                round(sparsevb_c_res[k,],3)))
      colnames(k_res)[1:2] <- c("Total_Beta", "Correct_Beta")
      rownames(k_res) <- c("UNHIDEM","LASSO","ALASSO","SCAD","MCP","SSLASSO", 
                           "VARBVS","SPARSEVB", "SPARSEVB_c")}
    
    cat("Iteration",k,"finished.\n")
    
    if(verbose){
      print(round(UNHIDEM_res[k,c(5:8)],3))
      if(ebreg_I){print(round(ebreg_res[k,8],3))}
      print(k_res)
    }
  }
  
  full_res <-   cbind(UNHIDEM_res, LASSO_res, adap_LASSO_res, SCAD_res, 
                      MCP_res, SSLASSO_res, varbvs_res, sparsevb_res, 
                      sparsevb_c_res, ebreg_res, M, M1, eta, sigma, sig_nois)
  
  write.csv(full_res,filname)
  
  
  MSE_comb <- full_res[,grepl('MSE', colnames(full_res))]
  MAD_comb <- full_res[,grepl('MAD', colnames(full_res))]
  beta_comb <- full_res[,grepl('Beta', colnames(full_res))]
  ECP_res <- full_res[,grepl('ECP', colnames(full_res))]
  cat("\n Summarizing some results:\n")
  
  cat("\n Train MSE, Test MSE and MAD of predictions of the true expectation,\n along with MSE of beta estimates:\n")
  avg_mse <- apply(MSE_comb,2,mean) ## 
  med_mad <- apply(MAD_comb,2,median) ### 
  avg_b_err <- apply(beta_comb,2,mean) 
  mse_mat <- matrix(avg_mse,nrow = length(avg_mse)/2, ncol = 2, byrow = TRUE)
  colnames(mse_mat) <- c("Train_MSE", "Test_MSE")
  if(ebreg_I){
    row.names(mse_mat) <- c("UNHIDEM","LASSO","ALASSO","SCAD","MCP","SSLASSO", 
                            "VARBVS","SPARSEVB", "SPARSEVB_c","EBREG")}
  if(!ebreg_I){
    rownames(mse_mat) <- c("UNHIDEM","LASSO","ALASSO","SCAD","MCP","SSLASSO", 
                           "VARBVS","SPARSEVB", "SPARSEVB_c")
  }
  pred_mat <- data.frame(mse_mat, MAD = med_mad, Beta_MSE = avg_b_err)
  print(pred_mat)
  
  cat("\n Average empirical coverage probabilities of 95% CI's and PI's:\n")
  print(apply(ECP_res,2,mean)) ### 
  
  full_res
}



sim_bin_LR_2 <- function(lat_data, err_sd, x, N, M, M1, sig_sp, lat_sp, eta_var, 
                         seed=NULL){
  
  if(is.null(seed)){seed = 238476}
  set.seed(seed)
  x1 <- x2 <- 1:as.integer(sqrt(M))
  signal_data <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp),x=x1,y=x2,grid=TRUE,n=1)
  
  quan <- quantile(array(signal_data),prob = 1-M1/M)
  signal_data <- 1*I(array(signal_data)>=quan)
  
  LP_data <- expand.grid(x1=x1,x2=x2)
  LP_data$signal     <- as.integer(signal_data)+1
  RE_eff <- rnorm(N)
  
  t_eta <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp,var = eta_var), x=x, y=x,grid=TRUE,n=1)
  t_data <- RFsimulate(model = RMstable(alpha = 2, scale = lat_sp,var = 1), x=x, y=x,grid=TRUE,n=N)
  
  RE_eff <- rnorm(N)
  t_datamat <- t(matrix(array(t_data),M,N))+RE_eff*err_sd
  Z_cont <- t(t(t_datamat)+c(array(lat_data)))
  Z <- 1*I(Z_cont < 0)
  
  signal<- seq(1,M,1)[LP_data$signal==2]
  sig_ind <- rep(0,M)
  sig_ind[signal] <- 1
  
  
  return(list(LP_data=LP_data,Z=Z, Z_cont = Z_cont,RE_eff=RE_eff,eta = t_eta,signal=signal,sig_ind=sig_ind))
}


adap_lasso <- function(Y,Z){
  ridge1_cv <- cv.glmnet(Z, Y, alpha = 0, nfolds = 10)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  cv.out <- cv.glmnet(Z, Y, alpha = 1, nfolds = 10, lambda.min.ratio = 0.001,
                      penalty.factor = 1 / abs(best_ridge_coef))
  
  while(cv.out$lambda.min==min(cv.out$lambda)){
    ### Further decreasing the lambda if the smallest lambda was the best.
    cv.out <- cv.glmnet(Z, Y, alpha = 1, nfolds = 10,
                        penalty.factor = 1 / abs(best_ridge_coef),
                        lambda=exp(c(log(cv.out$lambda.1se), 
                                     seq(log(min(cv.out$lambda)),
                                         log(min(cv.out$lambda)*0.001),
                                         length.out = 100))) )
  }
  cv.out
}

scad_func <- function(Y,Z,N,M){
  cv.out2 <- cv.ncvreg(matrix(as.numeric(Z),N,M),Y,family = "gaussian",penalty = "SCAD",lambda.min=0.01)
  while(cv.out2$lambda.min==min(cv.out2$lambda)){
    ### Further decreasing the lambda if the smallest lambda was the best.
    cv.out2 <- cv.ncvreg(matrix(as.numeric(Z),N,M),Y,family = "gaussian",penalty = "SCAD",lambda=exp(seq(log(min(cv.out2$lambda)),log(min(cv.out2$lambda)/100),length.out = 20)))
  }
  cv.out2
}

MCP_func <- function(Y,Z,N,M){
  cv.out3 <- cv.ncvreg(matrix(as.numeric(Z),N,M),Y,family = "gaussian",penalty = "MCP",lambda.min=0.01)
  while(cv.out3$lambda.min==min(cv.out3$lambda)){
    ### Further decreasing the lambda if the smallest lambda was the best.
    cv.out3 <- cv.ncvreg(matrix(as.numeric(Z),N,M),Y,family = "gaussian",penalty = "MCP",lambda=exp(seq(log(min(cv.out3$lambda)),log(min(cv.out3$lambda)/100),length.out = 20)))
  }
  cv.out3
}

sim_spat <- function(C, err_sd, beta0, N, M, M1, sig_sp, lat_sp, eta_var, 
                     seed=NULL){
  
  if(!is.null(seed)){set.seed(seed)}
  x1 <- x2 <- 1:as.integer(sqrt(M))
  signal_data <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp),x=x1,y=x2,grid=TRUE,n=1)
  quan <- quantile(array(signal_data),prob = 1-M1/M)
  signal_data <- 1*I(array(signal_data)>=quan)
  
  LP_data <- expand.grid(x1=x1,x2=x2)
  LP_data$signal     <- as.integer(signal_data)+1
  RE_eff <- rnorm(N)
  
  t_eta <- RFsimulate(model = RMstable(alpha = 2, scale = sig_sp,var = eta_var), x=x, y=x,grid=TRUE,n=1)
  t_data <- RFsimulate(model = RMstable(alpha = 2, scale = lat_sp,var = C^2), x=x, y=x,grid=TRUE,n=N)
  
  RE_eff <- rnorm(N)
  Z <- t(matrix(array(t_data),M,N))+RE_eff*err_sd
  
  signal<- seq(1,M,1)[LP_data$signal==2]
  sig_ind <- rep(0,M)
  sig_ind[signal] <- 1
  
  
  return(list(LP_data=LP_data,Z=Z,RE_eff=RE_eff,eta = t_eta,signal=signal,sig_ind=sig_ind))
}


