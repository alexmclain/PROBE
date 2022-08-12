#' @title Fitting PaRtitiOned empirical Bayes Ecm (PROBE) algorithm to sparse high-dimensional linear models.
#' @description A wrapper function for the main PROBE algorithm function.  The R package is a work in progress.
#' @usage probe(Y, Z, X, alpha, ep, maxit, Y_test, Z_test, X_test, verbose, 
#' signal, eta_i, plot_ind)
#' @param Y The outcome variable.
#' @param Z An \code{n x M} matrix of sparse predictors variables. 
#' @param X (optional) An \code{n x p} matrix or dataframe of other predictors not subjected to the sparsity assumption.
#' @param alpha Type I error; significance level
#' @param ep Value against which to compare convergence criterion (default = 0.1).
#' @param maxit Maximum number of iterations the algorithm will run for (default = 10000).
#' @param Y_test (optional) Test Y data used plotting purposes only (doesn't impact results).
#' @param Z_test (optional) Test Z data used plotting purposes only (doesn't impact results).
#' @param X_test (optional) Test X data used plotting purposes only (doesn't impact results).
#' @param verbose A logical (true/false) value whether to print algorithm iteration progress and summary quantities (default = FALSE).
#' @param signal (optional) A vector of indicies of the true non-null coefficients. This is used to calculate the true and false discovery rates by iteration for simulated data. Used plotting purposes only (doesn't impact results).
#' @param eta_i (optional) A vector of the true signal. This is used to calculate the MSE by iteration for simulated data. Used plotting purposes only (doesn't impact results).
#' @param plot_ind A logical values (True/False) for whether to output plots on algorithm results and progress (default = FALSE)
#' @return A list including 
#' 
#' \code{beta_ast_hat} MAP estimates of the regression coefficients (\eqn{\beta^\ast}),
#' 
#' \code{beta_hat, beta_hat_var} MAP estimates of the posterior expectation (beta_hat) and variance (beta_hat_var) of the prior mean (\eqn{\beta}) of the regression coefficients assuming \eqn{\gamma=1}, 
#' 
#' \code{gamma_hat} the posterior expectation of the latent \eqn{\gamma} variables, 
#' 
#' \code{sigma2_est} MAP estimate of the residual variance, 
#' 
#' \code{E_step} full results of the final E_step, 
#' 
#' \code{Calb_mod} results of first (\eqn{\alpha_0}) part of the M-step,  
#' 
#' \code{count} the total number of iterations before convergence. 
#' 
#' @seealso predict_probe_func to obtain predictions, credible intervals and prediction intervals from PROBE.
#' @examples
#' ### Example
#' data(Sim_data)
#' attach(Sim_data)
#' alpha <- 0.05
#' plot_ind <- TRUE
#' 
#' # Run the analysis. Y_test and Z_test are included for plotting purposes only
#' full_res <- probe( Y = Y, Z = Z, alpha = alpha, Y_test = Y_test, 
#' Z_test = Z_test, plot_ind = plot_ind)
#' 
#' # Predicting for test data
#' pred_res <- predict_probe_func(full_res, Z = Z_test, alpha = alpha)
#' head(pred_res)
#' 
#' # Estimate of the residual variance
#' full_res$sigma2_est
#' 
#' ### Example with additional covariate data X (not subjected to the sparsity assumption)
#' data(Sim_data_cov)
#' attach(Sim_data_cov)
#' 
#' # Calculating the true signal (the impact of Z only)
#' eta_i <- apply(t(Z)*beta_tr,2,sum) 
#' # Run the analysis. eta_i (true signal) and signal are included for plotting purposes only.
#' full_res <- probe( Y = Y, Z = Z, X = X, alpha = alpha, 
#' signal = signal, eta_i  = eta_i, plot_ind = plot_ind)
#' 
#' # Final estimates of the impact of X versus the true values:
#' data.frame(true_values = beta_X_tr, full_res$Calb_mod$res_data[-2,])
#' 
#' #Compare to a standard linear model of X on Y:
#' summary(lm(Y~X$Cont_cov + X$Binary_cov))$coefficients
#' 
#' 
probe <- function(Y, Z, X = NULL, alpha = 0.05, ep = 0.1, maxit = 10000, 
                    Y_test = NULL, Z_test = NULL, X_test = NULL,
                    verbose = FALSE, signal = NULL, eta_i = NULL, plot_ind = FALSE) {
  
  M <- dim(Z)[2]
  N <- dim(Z)[1]
  if(!is.null(X)){
    if (!is.matrix(X)) {
      X <- as.matrix(X)
    }
  }

  probe_func(Y = Y, Z = Z, X = X, alpha = alpha, verbose = verbose, 
                           signal = signal, maxit = maxit, eta_i = eta_i,
                           ep = ep, plot_ind = plot_ind, Y_test = Y_test, 
                           Z_test = Z_test, X_test = X_test)
}





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





mtr_func <- function(E_step, alpha, signal = NULL) {
  
  p_vals <- E_step$p_vals
  lfdr_val <- E_step$lfdr
  p_hat <- E_step$pi0
  M <- length(p_vals)
  alpha_hat <- alpha/p_hat
  
  T_R <- p.adjust(p_vals, method = "BY")
  R_BH <- p.adjust(p_vals, method = "BH")
  
  threshold <- 0
  lfdr_val[is.na(lfdr_val)] <- 1
  if (min(lfdr_val) < alpha) {
    threshold <- max(sort(lfdr_val)[cumsum(sort(lfdr_val)) < alpha])
  }
  
  BH_res <- data.frame(BY = 1 * I(T_R < alpha), BH = 1 * I(R_BH < alpha), 
                       LFDR = 1 * I(lfdr_val <= threshold))
  
  R_BY <- p.adjust(p_vals, method = "bonferroni")
  R2_BY <- p.adjust(p_vals, method = "holm")
  
  Bonf_res <- data.frame(Holm = 1 * I(R2_BY < alpha), Bonf = 1 * I(R_BY < alpha))
  
  BH_sum <- apply(BH_res,2,sum)
  Bonf_sum <- apply(Bonf_res,2,sum)
  if (!is.null(signal)) {
    n_signal <- !(1:M %in% signal)
    BH_sum <- c(sum(BH_res$BH), sum(BH_res$BH[signal]), sum(BH_res$BH[n_signal]))
    LFDR_res <- c(sum(BH_res$LFDR), sum(BH_res$LFDR[signal]), sum(BH_res$LFDR[n_signal]))
    BY_res <- c(sum(BH_res$BY), sum(BH_res$BY[signal]), sum(BH_res$BY[n_signal]))
    BH_sum <- data.frame(BY_sum = BY_res, LFDR_sum = LFDR_res, BH_sum = BH_sum)
    
    Bonf_sum <- c(sum(Bonf_res$Bonf), sum(Bonf_res$Bonf[signal]), sum(Bonf_res$Bonf[n_signal]))
    Holm_res2 <- c(sum(Bonf_res$Holm), sum(Bonf_res$Holm[signal]), 
                   sum(Bonf_res$Holm[n_signal]))
    
    Bonf_sum <- data.frame(Holm_sum = Holm_res2, Bonf_sum = Bonf_sum)
    row.names(BH_sum) <- c("Total", "Correct", "Errors")
    row.names(Bonf_sum) <- c("Total", "Correct", "Errors")
  }
  
  
  return(list(BH_res = BH_res, Bonf_res = Bonf_res, Bonf_sum = Bonf_sum, 
              BH_sum = BH_sum))
}


plot_probe_func <- function(full_res, test_plot, alpha) {
  
  plot_dat <- full_res$plot_dat
  
  a <- quantile(plot_dat$Pred_err, probs = 0.9)
  b <- min(plot_dat$Pred_err)
  L_pred <- plot_dat$Pred_err[length(plot_dat$Pred_err)]
  a <- max(c(a, b + (L_pred-b)*5/3 ))
  
  par(mar = c(4.2, 4.5, 0.5, 4.5), mfrow = c(1, 1))
  ylab_val <- expression(paste("Signal  ", MSE[t]))
  if (test_plot) {
    ylab_val <- expression(paste("Test  ", MSPE[t]))
  }
  
  plot(plot_dat$Iter, plot_dat$Pred_err, type = "l", ylim = c(b, a), 
       xlim = range(plot_dat$Iter), xlab = "Iteration", ylab = ylab_val, 
       lwd = 2, cex.lab = 1.4, cex.axis = 1.2, las = 1)
  
  p_vec <- plot_dat$Total_Disc
  b9 <- min(p_vec)
  a9 <- quantile(p_vec, probs = 0.9)
  
  axis(4, las = 1, at = seq(b, a, length.out = 4), 
       labels = round(seq(b9, a9, length.out = 4), 0), cex.axis = 1.2, las = 2)
  mtext("Number of rejections", side = 4, line = 2.9, cex = 1.2)
  trans_crit <- (p_vec - b9)/quantile(p_vec - b9, probs = 0.9) * (a - 
                                                                    b) + b
  col_vec <- rep(1, length(p_vec))
  if (!test_plot) {
    denom <- plot_dat$Total_Disc
    denom[denom == 0] <- 1
    col_vec <- ifelse(plot_dat$FP/denom < alpha, "grey60", 1)
  }
  points(plot_dat$Iter, trans_crit, col = col_vec, pch = 19, cex = 0.5)
  
  if (test_plot) { legend("topright", legend = c("Number of rejections \n using lfdr",
                                                 "Test MSPE"), lwd = 2, lty = c(0,1),
                          pch = c(20,-1), pt.lwd = c(2,0), cex = 1.3)
  }else{ legend("topright", legend = c("","Number of rejections \n using lfdr with:", 
                                       expression(FDR>alpha), expression(FDR<= alpha),
                                       "MSE of signal"),
                lwd = c(0,0,0,0,2), lty = c(1,1,3,3,1), col = c("white","white",1,"grey60",1),
                pch = c(-1,-1,20,20,-1), pt.lwd = c(0,0,2,2,0), cex = 1.3)}
}

report_func <- function(count, E_step, MTR_res, Xt_conv1, signal_track, 
                        report_pred) {
  
  CC_round <- signif(Xt_conv1, 2)
  if (is.null(signal_track)) {
    disc <- sum(MTR_res$BH_res$LFDR)
    if (!is.null(report_pred)) {
      report_pred <- signif(report_pred, 2)
      cat("Iteration=", count, "Number of discoveries (using lfdr)=", 
          disc, "Sum(delta)=", round(sum(E_step$delta), 1), " MSPE(test)=", 
          report_pred, "Convergence Crit=", CC_round, "\n")
    } else {
      cat("Iteration=", count, "Number of discoveries (using lfdr)=", 
          disc, "Sum(delta)=", round(sum(E_step$delta), 1), "Convergence Crit=", 
          CC_round, "\n")
    }
  } else {
    if (report_pred <= 1) {
      report_pred <- signif(report_pred, 3)
    }
    if (report_pred > 1) {
      report_pred <- round(report_pred, 1)
    }
    cat("Iteration=", count, " Hyp testing (using lfdr) TP=", signal_track[3], 
        " FP=", signal_track[1], " TN=", signal_track[2]-signal_track[4], " FN=", signal_track[4], 
        " MSE(signal)=", report_pred, " Convergence Crit=", CC_round, 
        "\n", sep = "")
  }
}

