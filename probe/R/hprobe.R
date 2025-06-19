#' @title Fitting PaRtitiOned empirical Bayes Ecm (PROBE) algorithm to sparse high-dimensional linear models with heterogeneous variance.
#' @description A wrapper function for the H-PROBE algorithm.
#' @usage hprobe(Y, X, V, Z = NULL, ep = 0.1, maxit = 10000, Y_test = NULL, X_test = NULL, 
#' Z_test = NULL, V_test = NULL, verbose = FALSE, signal = NULL, eta_i = NULL, alpha = 0.05, 
#' plot_ind = FALSE, adj = 5)
#' @param Y The outcome variable.
#' @param X An \code{n x M} matrix of sparse predictors variables. 
#' @param V A design matrix of predictors for the variance model (including an intercept).
#' @param Z (optional) An \code{n x p} matrix or dataframe of other predictors not subjected to the sparsity assumption.
#' @param ep Value against which to compare convergence criterion (default = 0.1).
#' @param maxit Maximum number of iterations the algorithm will run for (default = 10000).
#' @param Y_test (optional) Test Y data used plotting purposes only (doesn't impact results).
#' @param X_test (optional) Test X data used plotting purposes only (doesn't impact results).
#' @param Z_test (optional) Test Z data used plotting purposes only (doesn't impact results).
#' @param V_test (optional) Test V data used plotting purposes only (doesn't impact results).
#' @param verbose A logical (true/false) value whether to print algorithm iteration progress and summary quantities (default = FALSE).
#' @param signal (optional) A vector of indicies of the true non-null coefficients. This is used to calculate the true and false discovery rates by iteration for simulated data. Used plotting purposes only (doesn't impact results).
#' @param eta_i (optional) A vector of the true signal. This is used to calculate the MSE by iteration for simulated data. Used plotting purposes only (doesn't impact results).
#' @param alpha (optional) significance level
#' @param plot_ind A logical values (True/False) for whether to output plots on algorithm results and progress (default = FALSE)
#' @param adj Bandwidth parameter for empirical Bayes E-step. The bandwidth will be equal to \code{adj} times Silverman's 'rule of thumb' (default = 2).
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
#' @references \itemize{ \item McLain, AC, A Zgodic, H Bondell (2025). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. \textit{Computational Statistics and Data Analysis} 207, 108146.
#' \item Zgodic, A., Bai, R., Zhang, J., Wang, Y., Rorden, C., & McLain, A. (2023). Quantifying predictive uncertainty of aphasia severity in stroke patients with sparse heteroscedastic Bayesian high-dimensional regression. arXiv preprint arXiv:2309.08783.}
#' @examples
#' ### Example
#' data(h_Sim_data)
#' attach(h_sim_data)
#' 
#' # Run Analysis
#' res <- hprobe(Y = Y, X = X, V = V)
#'  
#' # Predicting for test data
#' pred_res <- predict_hprobe_func(res, X_test, V = V_test)
#' sqrt(mean((Y_test - pred_res$Pred)^2))
#' head(cbind(Y_test, pred_res))
#' 
#' plot(Y_test, pred_res$Pred, ylab = "Prediction", xlab = "Test Outcome")
#' abline(coef = c(0,1))
#' 
#' # Proportion of explained variance
#' 1 - var(Y_test - pred_res$Pred)/var(Y_test)
#' 
#' ## Omega coefficients (versus true values)
#' cbind(omega_tr, res$omega)
#' 
#' ## True versus estimated beta coeffiecients
#' plot(beta_tr, 
#'      res$beta_ast_hat, 
#'      xlab = "True Beta", 
#'      ylab = "Estimated Beta")
#' abline(coef = c(0,1))
#' 
#' ## Confusion matrix of true versus estimated signals using 0.5 cutoff.
#' table(beta_tr==0, res$gamma_hat<0.5)
#' 
#' 
#' @export
hprobe <- function(Y, X, V, Z = NULL, ep = 0.1, maxit = 10000, 
                   Y_test = NULL, X_test = NULL, Z_test = NULL,
                   V_test = NULL, #new
                   verbose = FALSE, signal = NULL, eta_i = NULL, 
                   alpha = 0.05, plot_ind = FALSE, adj = 5) {
  
  if(!is.null(Z)){
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
    }
  }
  
  hprobe_func(Y = Y, X = X, Z = Z, 
              V = V, #new
              alpha = alpha, verbose = verbose, 
              signal = signal, maxit = maxit, eta_i = eta_i,
              ep = ep, plot_ind = plot_ind, Y_test = Y_test, 
              X_test = X_test, Z_test = Z_test, 
              V_test = V_test, #new
              adj = adj)
}


hprobe_func <- function(Y, X, Z = NULL, 
                        V, #new
                        alpha, verbose = TRUE, signal, maxit = 1000,
                        eta_i = NULL, ep = 0.1, plot_ind = FALSE, 
                        Y_test = NULL, X_test = NULL, Z_test = NULL, 
                        V_test = NULL, #new
                        adj = 5){
  
  ##### Setting initial values and initializing outputs ####
  M <- dim(X)[2]
  N <- dim(X)[1]
  p <- 0
  if (!is.null(Z)) {
    p <- dim(Z)[2]
    beta_Z <- rep(0, p)
  }
  beta_t <- beta_var <- beta_tilde <- beta_tilde_var <- rep(0,M)
  gamma <- beta_t + 1
  W_ast <- rep(0, N)
  W_ast_var <- W_ast + 1
  
  # remove: sigma2 <- var(Y)
  sigma2 <- sigma2_O <- var(Y) #new
  Sigma_y_inv <- solve(diag(1, N)) #new
  
  count <- rcy_ct <- conv_check <- try2 <- 0
  plot_dat <- NULL
  X_2 <- X * X
  Xt_conv1 <- prev_Xt_conv1 <- 1
  T_vals <- signal_track <- report_pred <- plot_dat <- NULL
  
  while (count < maxit & conv_check < 1) {
    gamma_old <- gamma
    beta_t_old <- beta_t
    beta_var_old <- beta_var
    sigma2_old <- sigma2
    W_ast_old <- W_ast
    W_ast2_old <- W_ast_var
    count <- count + 1
    fact <- (count + 1)^(-1)
    
    # Performing the M-step.
    if (count == 1 & try2 == 0) {
      
      # remove: LR_update <- lr_cpp_func(Y, X, Z, sigma2)
      LR_update <- lr_cpp_func.h(Y, X, Z, sigma2, Sigma_y_inv) #new 
      
      beta_t_new <- c(LR_update$coef[,2])
      # Performing the damping step
      beta_t   <- beta_t*(1-fact) + beta_t_new*fact
      beta_var <- beta_var_old*(1-fact) + (LR_update$obs_SE[,2])^2*fact
    }else {
      # remove: LR_update <- m_step_cpp_func(Y, X, Z, W_ast, 
      #                              W_ast_var, gamma, 
      #                              beta_t, X_2, sigma2) 
      
      LR_update <- m_step_cpp_func.h(Y, X, Z, W_ast, W_ast_var, gamma,
                                     beta_t, X_2, sigma2, Sigma_y_inv) #new
      
      beta_t_new <- c(LR_update$coef[,2])
      # Performing the damping step
      beta_t   <- beta_t*(1-fact) + beta_t_new*fact
      beta_var <- 1/((1/beta_var_old)*(1-fact) + (1/LR_update$obs_SE[,2])^2*fact)
    }
    
    
    # Calculating the expectations of the gammas
    E_step <- e_step_func(beta_t, beta_var, df = N - 2, adj = adj, lambda = 0.1, monotone = TRUE ) 
    gamma <- E_step$gamma
    T_vals <- E_step$T_vals
    
    if (sum(gamma) == 0) {
      if(try2 == 0){
        warning("algorithm recycled back to null gamma values.\n 
            Trying again with different starting values. \n")
        rcy_ct <- count
        beta_t <- rep(1e-5,M)
        gamma <- rep(1,M)   
        T_vals <- NULL
        count <- try2 <- Xt_conv1 <- prev_Xt_conv1 <- 1 
      } else {
        E_step$beta_tilde <-  beta_t_old
        E_step$gamma <-  gamma_old
        E_step$beta_tilde_var <- beta_var_old
        break
      }
    }
    
    # Calculating the moments of W.
    W_W2_update <- m_update_func(X, X_2, beta_t, gamma)
    W_ast     <- W_W2_update$W_ast
    W_ast_var <- W_W2_update$W_ast_var
    if (var(W_ast) > 0) {
      # Run calibration model
      mod <- m_step_regression.h(Y, W_ast, W_ast_var + W_ast^2, Z, a = -3/2,
                                 Int = TRUE, V = V,
                                 Sigma_y_inv = Sigma_y_inv,
                                 c_param = 1000, sigma2_omega = Inf) #new
      
      Y_pred <- mod$Y_pred
      
      # remove: sigma2 <- mod$sigma2_est 
      omega <- mod$omega #new
      sigma2 <- min(c(mod$sigma2_est,sigma2_O)) #new
      
      Sigma_y <- diag(exp(-1*as.numeric(V%*%omega))) #new
      Sigma_y_inv <- inv_cpp(Sigma_y)$inverse #new
      t1 <- try(tmp <- solve(Sigma_y), silent = TRUE) #new
      if(unique(class(t1) %in% "try-error")){ print("Sigma_y_inv did not work") } #new #not necessary, was a check for me in simuls
      
      if (count > 1) {
        # Check Convergence
        if(length(c(W_ast_old[W_ast2_old>0],W_ast[W_ast2_old>0]))>0){
          Xt_conv1 <- pchisq(max((W_ast_old[W_ast2_old>0] - W_ast[W_ast2_old>0])^2 / 
                                   W_ast2_old[W_ast2_old>0])/log(N), df = 1)
        }
        if ( max(c(Xt_conv1,prev_Xt_conv1)) < ep) {conv_check <- conv_check + 1}
      }
    }
    prev_Xt_conv1 <- Xt_conv1
    
    # Getting prediction error if eta_i or test data given.
    report_pred <- mean((Y - Y_pred)^2)
    if (!is.null(eta_i)) {
      report_pred <- mean((eta_i - Y_pred)^2)
    }
    if (!is.null(Y_test) & !is.null(X_test)) {
      W_W2_test <- m_update_func(X_test, X_2 = NULL, beta_t, gamma)
      
      Y_pred_test <- cbind(1, W_W2_test$W_ast, Z_test) %*% mod$coef
      report_pred <- mean((Y_test - Y_pred_test)^2)
    }
    
    # Performing hypothesis testing on current estimates.
    MTR_res <- mtr_func(E_step, alpha, signal)
    if (!is.null(signal)) {
      signal_track <- c(MTR_res$BH_sum$LFDR_sum[3], 
                        M - length(signal) - MTR_res$BH_sum$LFDR_sum[3], 
                        MTR_res$BH_sum$LFDR_sum[2], 
                        length(signal) - MTR_res$BH_sum$LFDR_sum[2])
    }
    
    # Outputting results if verbose=TRUE
    if ((count %% 100) == 0 | conv_check == 1) {
      if(verbose){
        CC_round <- signif(Xt_conv1, 2)
        if (is.null(signal_track)) {
          disc <- sum(MTR_res$BH_res$LFDR)
          if (!is.null(report_pred)) {
            if (report_pred <= 1) {
              report_pred <- signif(report_pred, 3)
            }
            if (report_pred > 1) {
              report_pred <- round(report_pred, 1)
            }      
            message("Iteration=", count, "Number of discoveries (using lfdr)=", 
                    disc, "Sum(gamma)=", round(sum(E_step$gamma), 1), " MSPE(test)=", 
                    report_pred, "Convergence Crit=", CC_round, "\n")
          } else {
            message("Iteration=", count, "Number of discoveries (using lfdr)=", 
                    disc, "Sum(gamma)=", round(sum(E_step$gamma), 1), "Convergence Crit=", 
                    CC_round, "\n")
          }
        } else {
          if (report_pred <= 1) {
            report_pred <- signif(report_pred, 3)
          }
          if (report_pred > 1) {
            report_pred <- round(report_pred, 1)
          }
          message("Iteration=", count, " Hyp testing (using lfdr) TP=", signal_track[3], 
                  " FP=", signal_track[1], " TN=", signal_track[2]-signal_track[4], " FN=", signal_track[4], 
                  " MSE(signal)=", report_pred, " Convergence Crit=", CC_round, 
                  "\n", sep = "")
        }
      }
    }
    # Storing iteration data
    plot_dat <- rbind(plot_dat, c(rcy_ct+count, conv_check, signal_track, 
                                  sum(MTR_res$BH_res$LFDR), sum(gamma), 
                                  sigma2, report_pred, Xt_conv1))
  }
  
  beta_hat <- beta_t 
  beta_hat_var <- beta_var 
  beta_ast_hat <- beta_t*gamma*mod$coef[2]
  
  
  # Formatting iteration data for outputting
  if(!is.null(plot_dat)){
    if (!is.null(signal_track)) {
      colnames(plot_dat) <- c("Iter", "Conv_check", "FP", "TN", "TP", 
                              "FN", "Total_Disc", "Sum_gamma", "sigma2", "Pred_err", "Conv")
    } else {
      if (!is.null(report_pred)) {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc", 
                                "Sum_gamma", "sigma2", "Pred_err", "Conv")
      } else {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc", 
                                "Sum_gamma", "sigma2", "Conv")
      }
    }
    plot_dat <- data.frame(plot_dat)
  }
  
  conv <- 0
  if (conv_check == 0 & try2 < 2) {
    conv <- 1
    warning("convergence criteria not met. Set different convergence criteria or 
        raise maximum number of iterations.\n")
  }
  if (try2 == 2) {
    conv <- 2
    warning("algorithm recycled back to null values twice. Optimization failed.\n")
  }
  
  if(E_step$adj_warning != 0){
    message(paste0("Results indicate bandwidth in density estimation may be \n 
            too narrow, consider raising 'adj'. Current adj=",adj))
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
                   gamma_hat = gamma, E_step = E_step, Calb_mod = mod, count = count, plot_dat = plot_dat, 
                   Seq_test = Seq_test, M_step = M_step, sigma2_est = mod$sigma2_est, conv = conv, 
                   W_ast = W_ast, W_ast_var = W_ast_var, Y = Y, X = X, Z = Z, omega = omega, Sigma_y = Sigma_y)
  #new: in these outputs, I added "omega". "sigma2_est" can stay as is, but it will look different between H-PROBE and PROBE, as expected
  
  # Plotting iteration results if plot_ind=TRUE and either eta_i or test
  # data is given
  if (plot_ind & !is.null(plot_dat)) {
    if (is.null(report_pred)) {
      warning("cannot plot without true mean or test data.\n")
    } else {
      plot_probe_func(full_res, test_plot = !is.null(Y_test), alpha = alpha, signal = signal)
    }
  }
  
  
  return(full_res)
}

m_step_regression.h <- function(Y, W, W2, V, Sigma_y_inv, 
                                sigma2_omega = Inf, Z = NULL, 
                                a = -3/2, Int = TRUE, 
                                c_param = 1000) {
  
  N <- length(Y)
  if(!is.null(W)){
    if(Int){ Wmat <- cbind(1, W, Z) }
    if(!Int & !is.null(Z) ){ Wmat <- cbind(W, Z) }
    if(!Int & is.null(Z) ){ Wmat <- matrix(W) }
    W_col <- c(1 + 1*I(Int))
    
    # remove: WWpr  <- t(Wmat) %*% Wmat
    # remove WWpr[W_col, W_col] <- sum(W2)
    
    WWpr  <- t(Wmat) %*% Sigma_y_inv %*% Wmat #new
    WWpr[W_col, W_col] <- sum(W2 %*% Sigma_y_inv) #new
    
    df <- N - ncol(Wmat) + 2*a - 3
    
    if (det(WWpr) != 0) {
      WWpr_inv <- solve(WWpr)
      # remove: beta_w <- WWpr_inv %*% t(Wmat) %*% Y
      beta_w <- WWpr_inv %*% t(Wmat) %*% Sigma_y_inv %*% Y #new
      
      Y_pred <- Wmat %*% beta_w
      hat <- diag(Wmat %*% WWpr_inv %*% t(Wmat))
      resid <- (Y - Y_pred)
      RSS <- sum(resid^2) + beta_w[W_col]^2*(sum(W2) - sum(W^2))
      
      # remove: VCV <- RSS/(df) * t(WWpr_inv) %*% (t(Wmat) %*% Wmat) %*% WWpr_inv
      VCV <- t(WWpr_inv) %*% (t(Wmat) %*% Sigma_y_inv %*% Wmat) %*% WWpr_inv #new
    }
    if (det(WWpr) == 0) {
      final_mod <- lm(Y ~ 0 + Wmat)
      beta_w <- as.numeric(final_mod$coefficients)
      beta_w[is.na(beta_w)] <- 0
      Y_pred <- as.numeric(final_mod$fitted.values)
      hat <- as.numeric(influence(final_mod)$hat)
      resid <- as.numeric(final_mod$residuals)
      RSS <- sum(resid^2)
      VCV <- vcov(final_mod)
    }
    
    Std_Err  <- sqrt(diag(VCV))
    T_vals   <- beta_w/Std_Err
  }
  if(is.null(W)){
    if(is.null(Z)){
      fail_lm <- lm(Y~1)
      coef_fail <- summary(fail_lm)$coefficients
      Y_pred <- as.numeric(fail_lm$fitted.values)
      hat <- as.numeric(influence(fail_lm)$hat)
      resid <- as.numeric(fail_lm$residuals)
      RSS <- sum(resid^2)
      df <- N- 1 - + 2*a - 3
      beta_w <- c(coef_fail[1], 0)
      VCV <- matrix(0,2,2)
      VCV[1,1] <- vcov(fail_lm)
      Std_Err  <- c(sqrt(vcov(fail_lm)), 0)
      T_vals   <- c(beta_w[1]/Std_Err[1], 0)
    }else{
      Wmat <- cbind(1, Z)
      fail_lm <- lm(Y~Wmat)
      coef_fail <- summary(fail_lm)$coefficients
      Y_pred <- as.numeric(fail_lm$fitted.values)
      hat <- as.numeric(influence(fail_lm)$hat)
      resid <- as.numeric(fail_lm$residuals)
      RSS <- sum(resid^2)
      p <- dim(Wmat)
      df <- N - p + 2*a - 3
      beta_w <- c(coef_fail[1], 0, coef_fail[-1])
      VCV <- matrix(0,2+p,2+p)
      VCV[-2,-2] <- vcov(fail_lm)
      Std_Err  <- sqrt(diag(VCV))
      Std_Err  <- c(Std_Err[1],0,Std_Err[-1])
      T_vals   <- c(beta_w[1]/Std_Err[1], 0, beta_w[-1]/Std_Err[-1])
    }
  }
  
  p_val    <- pt(abs(T_vals), df = df, lower.tail = FALSE) * 2
  res_data <- data.frame(Estimate = beta_w, Std_Err = Std_Err, T_val = T_vals, df = df, p_val = p_val)
  if(Int){row.names(res_data) <- c("Intercept", "W", colnames(Z))}
  if(!Int){row.names(res_data) <- c("W", colnames(Z))}
  
  # new: this big chunk of comment
  # optimization for the posterior of omega
  log.lklh.MLG <- function(par, V, Y, alphaW, Wmat, c_param, sigma2_omega){
    (-1) * ( sum( 0.5*V%*%par - 0.5*((Y-Wmat%*%alphaW)^2)*exp(V%*%par) ) + ( (c_param*rep(1, ncol(V))*(c_param^(-0.5))*(1/sigma2_omega)) %*% (diag(ncol(V))%*%par) ) -
               ( c_param*rep(1, ncol(V)) %*% exp((c_param^(-0.5))*(1/sigma2_omega)*diag(ncol(V))%*%par) ) )
  }
  
  optim_results <- optim(par = rep(0,ncol(V)) , fn = log.lklh.MLG,
                         V = V, Y = Y, alphaW = beta_w,
                         Wmat = Wmat, c_param = c_param,
                         sigma2_omega = sigma2_omega, method = "BFGS")
  
  omega <- optim_results$par
  
  # new: added 'omega = omega' to this list of outputs
  list(coef = beta_w, sigma2_est = RSS/df, RSS = RSS, 
       Y_pred = Y_pred, hat = hat, resid = resid, 
       VCV = VCV, res_data = res_data, omega = omega)
  
}

# new: added a new argument for this function 'Sigma_y_inv' 
m_step_cpp_func.h <- function(Y, X, Z = NULL, W_ast, W_ast_var, gamma, 
                              beta_vec, X_2, sigma2, Sigma_y_inv) {
  
  N <- length(Y)
  if (!is.null(Z)) {
    # remove: LRcpp <- PROBE_cpp0_5_6_covs(Y, X, W_ast, W_ast_var, gamma, 
    #                              beta_vec, X_2, sigma2, as.matrix(Z))
    LRcpp <- PROBE_cpp0_5_6_covs_h(Y, X, W_ast, W_ast_var, gamma,
                                   beta_vec, X_2, sigma2, as.matrix(Z), Sigma_y_inv) #new
  } else {
    # remove: LRcpp <- PROBE_cpp0_5_6(Y, X, W_ast, W_ast_var, gamma, 
    #                         beta_vec, X_2, sigma2)
    LRcpp <- PROBE_cpp0_5_6_h(Y, X, W_ast, W_ast_var, gamma,
                              beta_vec, X_2, sigma2, Sigma_y_inv) #new
  }
  
  ret <- list(coef = LRcpp$Coefficients, obs_SE = LRcpp$StdErr)
  
  return(ret)
}

# new: a new argument 'Sigma_y_inv' was added to this function
lr_cpp_func.h <- function(Y, X, Z = NULL, sigma2, Sigma_y_inv) {
  
  N <- length(Y)
  if (!is.null(Z)) {
    # remove: LRcpp <- LM_w_COVS_by_col(Y, X, as.matrix(Z), sigma2)
    LRcpp <- LM_w_COVS_by_col_h(Y, X, as.matrix(Z), sigma2, Sigma_y_inv) #new
  } else {
    # remove: LRcpp <- LM_by_col(Y, X, sigma2)
    LRcpp <- LM_by_col_h(Y, X, sigma2, Sigma_y_inv) #new
  }
  t_val <- LRcpp$Coefficients[, 2]/LRcpp$StdErr[, 2]
  
  ret <- list(coef = LRcpp$Coefficients, obs_SE = LRcpp$StdErr, T_val = t_val)
  return(ret)
}


m_update_func <- function(X,X_2,beta_tilde, gamma, beta_tilde_var=0){
  X_gamma <- MVM(X,gamma*beta_tilde)$Res 
  W_ast<- c(Row_sum(as.matrix(X_gamma))$Rowsum)
  W_ast[is.nan(W_ast) | is.na(W_ast)] <- 0
  W_ast_var <- NULL
  if(!is.null(X_2)){
    X_gamma2 <- MVM(X_2,beta_tilde^2*gamma*(1-gamma) + gamma*beta_tilde_var)$Res 
    W_ast_var <- c(Row_sum(as.matrix(X_gamma2))$Rowsum)
    W_ast_var[is.nan(W_ast_var) | is.na(W_ast_var)] <- 0
  }
  return(list(W_ast=W_ast,W_ast_var=W_ast_var))
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


plot_probe_func <- function(full_res, test_plot, alpha, signal) {
  
  plot_dat <- full_res$plot_dat
  
  a <- quantile(plot_dat$Pred_err, probs = 0.9)
  b <- min(plot_dat$Pred_err)
  L_pred <- plot_dat$Pred_err[length(plot_dat$Pred_err)]
  a <- max(c(a, b + (L_pred-b)*5/3 ))
  
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
  
  if (test_plot) { 
    points(plot_dat$Iter, trans_crit, col = 1, pch = 19, cex = 0.5)
    legend("topright", legend = c("Rejections with lfdr",
                                  "Test MSPE"), lwd = 2, lty = c(0,1),
           pch = c(19,-1), pt.lwd = c(2,0), cex = 1.3)
  }else{ 
    if(!is.null(signal)){
      denom <- plot_dat$Total_Disc
      denom[denom == 0] <- 1
      col_vec <- ifelse(plot_dat$FP/denom < alpha, "grey60", 1)
      points(plot_dat$Iter, trans_crit, col = col_vec, pch = 19, cex = 0.5)
      legend("topright", legend = c(expression(FDR>alpha), expression(FDR<= alpha),
                                    "MSE of signal"),
             lwd = c(0,0,2), lty = c(3,3,1), col = c(1,"grey60",1),
             pch = c(20,20,-1), pt.lwd = c(2,2,0), cex = 1.3)
    }else{
      points(plot_dat$Iter, trans_crit, col = 1, pch = 19, cex = 0.5)
      legend("topright", legend = c("Rejections with lfdr",
                                    "MSE of signal"),
             lwd = c(2), lty = c(0,1), col = 1,
             pch = c(19,-1), pt.lwd = c(2,0), cex = 1.3)
    }
  }
  
}



