probe_func <- function(Y, Z, X = NULL, alpha, verbose = TRUE, signal, maxit = 1000, 
                         eta_i = NULL, ep = 0.1, plot_ind = FALSE, 
                         Y_test = NULL, Z_test = NULL, X_test = NULL){
  
  ##### Setting initial values and initializing outputs ####
  M <- dim(Z)[2]
  N <- dim(Z)[1]
  p <- 0
  if (!is.null(X)) {
    p <- dim(X)[2]
    beta_X <- rep(0, p)
  }
  beta_t <- beta_var <- beta_tilde <- beta_tilde_var <- rep(0,M)
  T_vals <- NULL
  W_ast <- rep(0, N)
  delta <- beta_t + 1
  df <- N
  sigma2 <- sigma2_O <- var(Y)
  
  W_ast_var <- W_ast + 1
  count <- 0
  conv_check <- 0
  CC_count <- 0
  plot_dat <- NULL
  Z_2 <- Z * Z
  Xt_conv1 <- Xt_conv2 <- Xt_conv3 <- 1
  try2 <- 0
  signal_track <- report_pred <- NULL
  tau_s <- 0
  rcy_ct <- count
  
  while (count < maxit & conv_check < 1) {
    delta_old <- delta
    beta_t_old <- beta_t
    beta_var_old <- beta_var
    sigma2_old <- sigma2
    tau_s_old <- tau_s
    W_ast_old <- W_ast
    W_ast2_old <- W_ast_var
    count <- count + 1
    fact <- (count + 1)^(-1)
    
    # Performing the M-step.
    if (count == 1 & try2 == 0) {
      LR_update <- lr_cpp_func(Y, Z, X, sigma2)
      beta_t_new <- c(LR_update$coef[,2])
      # Performing part (a) of the E-step.
      beta_t   <- beta_t*(1-fact) + beta_t_new*fact
      beta_var <- beta_var_old*(1-fact) + (LR_update$obs_SE[,2])^2*fact
    }else {
      LR_update <- m_step_cpp_func(Y, Z, X, W_ast, 
                                   W_ast_var, delta, 
                                   beta_t, Z_2, sigma2) 
      
      beta_t_new <- c(LR_update$coef[,2])
      # Performing part (a) of the E-step.
      beta_t   <- beta_t*(1-fact) + beta_t_new*fact
      beta_var <- 1/((1/beta_var_old)*(1-fact) + (1/LR_update$obs_SE[,2])^2*fact)
    }
    
    
    
    # Performing part (b) of the E-step.
    E_step <- e_step_func(beta_t, beta_var, df = N - 2, adj = 5, lambda = 0.1, monotone = TRUE ) 
    delta <- E_step$delta
    T_vals <- E_step$T_vals
    
    if (sum(delta) == 0) {
      E_step <- e_step_func(beta_t, beta_var, df = N - 2, adj = 5, lambda = 0.1) 
      delta <- E_step$delta
      T_vals <- E_step$T_vals
      
      if (sum(delta) == 0) {
        if(try2 == 0){
          cat("Warning loop completely recycled back to beta=0.\n 
            Trying again with different starting values. \n")
          rcy_ct <- count
          count <- 0
          beta_t <- rep(1e-5,M)
          beta_var <- rep(0,M)
          delta <- rep(1,M)   
          T_vals <- NULL
          Xt_conv1 <- 1 
          try2 <- 1
          tau_s <- 1/sigma2
        } else {
          cat("Warning loop completely recycled back to beta=0 again.\n")
          break
        }
      }
    }
    
    
    
    # Performing part (c) of the E-step.
    W_W2_update <- m_update_func(Z, Z_2, beta_t, delta)
    W_ast     <- W_W2_update$W_ast
    W_ast_var <- W_W2_update$W_ast_var
    if (var(W_ast) > 0) {
      # Run calibration model
      mod <- m_step_regression(Y, W_ast, W_ast_var + W_ast^2, X, a = -3/2, Int = TRUE) 
      Y_pred <- mod$Y_pred
      sigma2 <- min(c(mod$sigma2_est,sigma2_O)) 
      if (count > 1) {
        # Check Convergence
        if(length(c(W_ast_old[W_ast2_old>0],W_ast[W_ast2_old>0]))>0){
          Xt_conv1 <- pchisq(max((W_ast_old[W_ast2_old>0] - W_ast[W_ast2_old>0])^2 / 
                                   W_ast2_old[W_ast2_old>0])/log(N), df = 1)
        }
        if ( Xt_conv1 < ep) {conv_check <- conv_check + 1}
      }
    }
    
    # Getting prediction error if eta_i or test data given.
    if (!is.null(eta_i)) {
      report_pred <- mean((eta_i - Y_pred)^2)
    }
    if (!is.null(Y_test) & !is.null(Z_test)) {
      W_W2_test <- m_update_func(Z_test, Z_2 = NULL, beta_t, delta)
      
      Y_pred_test <- cbind(1, W_W2_test$W_ast, X_test) %*% mod$coef
      report_pred <- mean((Y_test - Y_pred_test)^2)
    }
    
    # Performing hypothesis testing on current estimates.
    MTR_res <- mtr_func(E_step, alpha, signal)
    if (!is.null(signal)) {
      signal_track <- c(MTR_res$BH_sum$LFDR_sum[3], 
                        M - MTR_res$BH_sum$LFDR_sum[1], 
                        MTR_res$BH_sum$LFDR_sum[2], 
                        length(signal) - MTR_res$BH_sum$LFDR_sum[2])
    }
    
    # Outputting results if verbose=TRUE
    if (( (count %% 100) == 0 | conv_check == 1) & verbose) {
      if(count !=0){
        report_func(count, E_step, MTR_res, Xt_conv1, signal_track, 
                    report_pred)
      }
    }
    # Storing iteration data
    plot_dat <- rbind(plot_dat, c(rcy_ct+count, conv_check, signal_track, 
                                  sum(MTR_res$BH_res$LFDR), sum(delta), 
                                  sigma2, report_pred, Xt_conv1))
  }
  
  beta_hat <- beta_t 
  beta_hat_var <- beta_var 
  beta_ast_hat <- beta_t*delta*mod$coef[2]
  
  
  # Formatting iteration data for outputting
  if(!is.null(plot_dat)){
    if (!is.null(signal_track)) {
      colnames(plot_dat) <- c("Iter", "Conv_check", "FP", "TN", "TP", 
                              "FN", "Total_Disc", "Sum_delta", "sigma2", "Pred_err", "Conv")
      
    } else {
      if (!is.null(report_pred)) {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc", 
                                "Sum_delta", "sigma2", "Pred_err", "Conv")
      } else {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc", 
                                "Sum_delta", "sigma2", "Conv")
      }
    }
    plot_dat <- data.frame(plot_dat)
  }
  
  conv <- 0
  if (conv_check == 0 & try2 < 2) {
    conv <- 1
    cat("Warning: convergence criteria not met. Set different convergence criteria or 
        raise maximum number of iterations.\n")
  }
  if (try2 == 2) {
    conv <- 2
    cat("Warning: loop completely recycled back to beta=delta=0 twice. Optimization failed.\n")
  }
  
  M_step <- LR_update
  Seq_test <- NULL
  if(length(beta_t_new[beta_tilde_var>0]) > 0){
    Seq_test <- pchisq(sum((beta_t_new[beta_tilde_var>0] - beta_tilde[beta_tilde_var>0])^2 / 
                             beta_tilde_var[beta_tilde_var>0]),
                       df=length(beta_t_new[beta_tilde_var>0]),
                       lower.tail = FALSE)
  }
  
  full_res <- list(beta_ast_hat = beta_ast_hat, beta_hat = beta_hat, beta_hat_var = beta_hat_var, 
                   gamma_hat = delta, E_step = E_step, Calb_mod = mod, count = count, plot_dat = plot_dat, 
                   Seq_test = Seq_test, M_step = M_step, sigma2_est = mod$sigma2_est, conv = conv, 
                   W_ast = W_ast, W_ast_var = W_ast_var, Y = Y, Z = Z, X = X)
  
  # Plotting iteration results if plot_ind=TRUE and either eta_i or test
  # data is given
  if (plot_ind & !is.null(plot_dat)) {
    if (is.null(report_pred)) {
      cat("Warning: cannot plot without eta_i or test data.\n")
    } else {
      plot_probe_func(full_res, test_plot = !is.null(Y_test), alpha = alpha)
    }
  }
  
  if (conv_check == 0) {
    cat("Warning: convergence criteria not met. Set different convergence criteria or 
        raise maximum number of iterations.\n")
  }
  
  
  return(full_res)
}