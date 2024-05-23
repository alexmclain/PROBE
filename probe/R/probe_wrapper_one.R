#' @title Fitting PaRtitiOned empirical Bayes Ecm (PROBE) algorithm to sparse high-dimensional linear models.
#' @description A wrapper function for the one-at-a-time variant of the PROBE algorithm.
#' @usage probe_one(Y, X, ep = 0.001, maxit = 10000, Y_test = NULL, X_test = NULL, 
#' verbose = FALSE, signal = NULL, eta_i = NULL, alpha = 0.05, plot_ind = FALSE, 
#' order.method = "lasso", adj = 10, delta = 0.4, update_order= NULL, beta_start= NULL)
#' 

#' @param Y The outcome variable.
#' @param X An \code{n x M} matrix of sparse predictors variables. 
#' @param ep Value against which to compare convergence criterion (default = 0.001).
#' @param maxit Maximum number of iterations the algorithm will run for (default = 10000).
#' @param Y_test (optional) Test Y data used plotting purposes only (doesn't impact results).
#' @param X_test (optional) Test X data used plotting purposes only (doesn't impact results).
#' @param verbose A logical (true/false) value whether to print algorithm iteration progress and summary quantities (default = FALSE).
#' @param signal (optional) A vector of indicies of the true non-null coefficients. This is used to calculate the true and false discovery rates by iteration for simulated data. Used plotting purposes only (doesn't impact results).
#' @param eta_i (optional) A vector of the true signal. This is used to calculate the MSE by iteration for simulated data. Used plotting purposes only (doesn't impact results).
#' @param alpha (optional) significance level
#' @param plot_ind A logical values (True/False) for whether to output plots on algorithm results and progress (default = FALSE)
#' @param order.method Updating order and initial values of the algorithm. For \code{lasso} (default) or \code{ridge}, a lasso or a ridge regression model (fit with 10-fold CV) will be fitted and used. The \code{update_order} is defined by the absolute values of the coefficient and \code{beta_start} is the coefficient values. When using \code{none}, \code{update_order} and \code{beta_start} must be given. \code{random} will randomly select the updating order and use very small values for \code{beta_start}. 
#' @param adj Bandwidth parameter for empirical Bayes E-step. The bandwidth will be equal to \code{adj} times Silverman's 'rule of thumb' (default = 10).
#' @param delta Learning rate for iteration t is (1 + t)^(-1 + delta) (default delta = 0.4).
#' @param update_order Manual value for the updating order for when \code{order.method = "none"} is used.
#' @param beta_start Manual value for the starting beta coefficients for when \code{order.method = "none"} is used.  
#' @param seed Seed value to ensure reproducibility when \code{order.method = "lasso"}, \code{order.method = "ridge"}, or \code{order.method = "random"}. 
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
#' \code{count} the total number of iterations before convergence. 
#' 
#' @seealso predict_probe_func to obtain predictions.
#' @references 
#' McLain, A. C., Zgodic, A., & Bondell, H. (2022). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2209.08139..
#' @examples
#' ### Example
#' data(Sim_data)
#' data(Sim_data_test)
#' attach(Sim_data)
#' attach(Sim_data_test)
#' plot_ind <- TRUE
#' adj <- 10
#' 
#' # Run the analysis. Y_test and X_test are included for plotting purposes only
#' full_res <- probe_one( Y = Y, X = X, Y_test = Y_test, order.method = "lasso",
#' X_test = X_test, plot_ind = plot_ind, adj = adj)
#' 
#' # Predicting for test data
#' pred_res <- predict_probe_func(full_res, X = X_test)
#' sqrt(mean((Y_test - pred_res$Pred)^2))
#' 
#' # Estimate of the residual variance and true value
#' full_res$sigma2_est
#' sigma2_tr
#' 
#' # RMSE of estimated beta coefficients
#' beta_ast_est <- c(full_res$beta_ast_hat)
#' sqrt(mean((beta_ast_est - beta_tr)^2))
#' 
#' # Posterior expectation of gamma by true
#' gamma_est <- full_res$E_step$gamma
#' table(gamma_est > 0.5, beta_tr > 0)
#' sum(gamma_est)
#' sum(gamma_est[beta_tr>0])
#' 
#' 
#' @export
probe_one <- function(Y, X, ep = 0.001, maxit = 10000, Y_test = NULL, X_test = NULL, 
                      verbose = FALSE, signal = NULL, eta_i = NULL, alpha = 0.05, plot_ind = FALSE, 
                      order.method = "lasso", adj = 10, delta = 0.4, update_order= NULL, 
                      beta_start= NULL) {
  
  ### The one at a time probe.
  M <- dim(X)[2]
  N <- dim(X)[1]
  
  ### Centering the data.
  X_0 <- X
  Y_0 <- Y
  Y_mn <- 0
  X_mean <- 0
  X_mean <- apply(X,2,mean)
  X <- t(t(X) - X_mean)
  Y_mn <- mean(Y)
  Y <- Y - Y_mn
  if(!is.null(eta_i)){eta_i <- eta_i - Y_mn}
  if(!is.null(Y_test)){Y_test <- Y_test - mean(Y_test)}
  if(!is.null(X_test)){
    X_test <- t(t(X_test) - X_mean)
  }
  
  if(order.method == "none" & 
     (is.null(update_order) | is.null(beta_start)) ){
    stop("update_order and beta_start must be given for update.order='none'.")
  }
  if(order.method == "lasso"){
    set.seed(546646+floor(max(Y^2)))
    foldid = sample(rep(1:10,ceiling(N/10)),N)
    r.out <- cv.glmnet(X,Y,alpha = 1,lambda.min.ratio=0.001, foldid = foldid)
    beta_start <- coef(r.out,s="lambda.min")[-1]
    update_order = order(abs(coef(r.out,s="lambda.min")[-1]), decreasing = TRUE)
  }
  if(order.method == "ridge"){
    set.seed(546646+floor(max(Y^2)))
    foldid = sample(rep(1:10,ceiling(N/10)),N)
    r.out <- cv.glmnet(X,Y,alpha = 0,lambda.min.ratio=0.001, foldid = foldid)
    beta_start <- coef(r.out,s="lambda.min")[-1]
    update_order = order(abs(coef(r.out,s="lambda.min")[-1]), decreasing = TRUE)
  }
  if(order.method == "random"){
    set.seed(546646+floor(max(Y^2)))
    update_order <- sample(1:M,M)
    beta_start <- rep(0.0001,M)
    
  }
  probe_func_one(Y = Y, X = X, alpha = alpha, verbose = verbose, 
                 signal = signal, maxit = maxit, eta_i = eta_i,
                 ep = ep, plot_ind = plot_ind, Y_mn = Y_mn, X_mean = X_mean, 
                 Y_test = Y_test, X_test = X_test, update_order = update_order, 
                 beta_start = beta_start, adj = adj, delta = delta)
}


probe_func_one <- function(Y, X, alpha, verbose = TRUE, signal, maxit = 1000, 
                           eta_i = NULL, ep = 0.1, plot_ind = FALSE, 
                           Y_mn = 0, X_mean = NULL, Y_test = NULL, X_test = NULL, 
                           update_order = NULL, adj = 10, beta_start, delta = 0.5){
  
  ##### Setting initial values and initializing outputs ####
  M <- dim(X)[2]
  N <- dim(X)[1]
  beta_t <- beta_tilde <- beta_start
  beta_var <- beta_tilde_var <- rep(0,M)
  T_vals <- NULL
  W_ast <- rep(0, N)
  gamma <- rep(1,M)
  df <- N
  sigma2 <- v_Y <- var(Y)
  
  W_ast_var <- W_ast + 1
  count <- 0
  conv_check <- 0
  CC_count <- 0
  plot_dat <- NULL
  X_2 <- X * X
  Xt_conv1 <- prev_Xt_conv1 <- 1
  try2 <- 0
  signal_track <- report_pred <- NULL
  tau_s <- 0
  rcy_ct <- count
  alpha_est <- 1
  
  W_W2_update <- m_update_func(X, X_2, beta_t, gamma)
  W_ast     <- W_W2_update$W_ast
  W_ast_var <- W_W2_update$W_ast_var
  while (count < maxit & conv_check < 1) {
    alpha_est_old <- alpha_est
    gamma_old <- gamma
    beta_t_old <- beta_t
    beta_var_old <- beta_var
    sigma2_old <- sigma2
    tau_s_old <- tau_s
    W_ast_old <- W_ast
    W_ast2_old <- W_ast_var
    count <- count + 1
    
    LR_update <- m_step_cpp_one(Y, X, W_ast, W_ast_var, 
                                  gamma, beta_t, X_2, sigma2, 
                                  update_order-1)
    
    sigma2 <- LR_update$Sigma2
    beta_t_new   <- LR_update$Coefficients[,1]
    obs_var <- LR_update$StdErr^2
    fact <- (count + 1)^(-1 + delta)
    
    beta_t   <- beta_t*(1-fact) + beta_t_new*fact
    if(count > 1){
      beta_var <- 1/((1/beta_var_old)*(1-fact) + (1/obs_var)*fact)
    }else{
      beta_var <- (obs_var)*fact
    }
    
    alpha_est <- LR_update$Coefficients_alpha[M+1]
    
    E_step <- e_step_func(beta_t, beta_var, df = N - 2, adj = adj, 
                          lambda = 0.1, monotone = TRUE )
    
    gamma <- E_step$gamma
    T_vals <- E_step$T_vals
    
    
    if (sum(gamma) == 0) {
      E_step <- e_step_func(beta_t, beta_var, df = N - 2, adj = 5, lambda = 0.1) 
      gamma <- E_step$gamma
      T_vals <- E_step$T_vals
      
      if (sum(gamma) == 0) {
        if(try2 == 0){
          cat("Warning loop completely recycled back to beta=0.\n 
            Trying again with different starting values. \n")
          rcy_ct <- count
          count <- 0
          beta_t <- rep(1e-5,M)
          beta_var_old <- beta_var_old*0
          gamma <- rep(1,M)   
          T_vals <- NULL
          Xt_conv1 <- prev_Xt_conv1 <- 1 
          try2 <- 1
          tau_s <- 1/sigma2
        } else {
          cat("Warning loop completely recycled back to beta=0 again.\n")
          break
        }
      }
    }
    
    
    # Updating expectations of W.
    W_W2_update <- m_update_func(X, X_2, beta_t, gamma)
    W_ast     <- W_W2_update$W_ast
    W_ast_var <- W_W2_update$W_ast_var
    if (var(W_ast) > 0) {
      # Remapping all parameters
      mod <- m_step_regression(Y, W_ast, W_ast_var + W_ast^2, a = -3/2, Int = FALSE)
      alpha_est <- mod$coef[1]
      beta_t <- alpha_est*beta_t
      beta_var <- beta_var*alpha_est^2
      W_ast     <- W_ast*alpha_est
      W_ast_var <- W_ast_var*(alpha_est^2)
      
      Y_pred <- mod$Y_pred
      if (count > 1) {
        # Check Convergence
        if(length(c(W_ast_old[W_ast2_old>0],W_ast[W_ast2_old>0]))>0 &
           (max(W_ast_var)>0) ){
          Xt_conv1 <- pchisq(max((W_ast_old[W_ast2_old>0] - W_ast[W_ast2_old>0])^2 / 
                                   W_ast2_old[W_ast2_old>0])/log(N), df = 1)
        }else{
          Xt_conv1 = max(abs(c(beta_t-beta_t_old, gamma - gamma_old)))*10
        }
        if(sum(gamma[-which.max(gamma)]) == 0){
          Xt_conv1 <- 0
          try2 <- 3
        }
        if ( max(c(Xt_conv1,prev_Xt_conv1)) < ep) {conv_check <- conv_check + 1}
      }
    }
    prev_Xt_conv1 <- Xt_conv1
    
    # Getting prediction error if eta_i or test data given.
    report_pred <- mean((Y - W_ast)^2)
    if (!is.null(eta_i)) {
      report_pred <- mean((eta_i - W_ast)^2)
    }else{
      if (!is.null(Y_test) & !is.null(X_test)) {
        W_W2_test <- m_update_func(X_test, X_2 = NULL, beta_t, gamma)
        Y_pred_test <- W_W2_test$W_ast
        report_pred <- mean((Y_test - Y_pred_test)^2)
      }
    }
    
    # Performing hypothesis testing on current estimates.
    MTR_res <- mtr_func(E_step, alpha , signal)
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
                                  sigma2, report_pred, Xt_conv1,alpha_est))
  }
  gamma <- c(gamma)
  beta_hat <- beta_t
  beta_hat_var <- beta_var 
  beta_ast_hat <- c(beta_t*gamma)
  mod <- m_step_regression(Y, W_ast, W_ast_var + W_ast^2, a = -3/2, Int = TRUE)
  mod$coef[1] <- Y_mn
  mod$Y_pred <- mod$Y_pred + Y_mn
  mod$VCV[1,1] <- v_Y/N
  mod$VCV[1,2] <- mod$VCV[2,1] <- 0
  mod$res_data[1,] <- c(Y_mn, sqrt(mod$VCV[1,1]), Y_mn/sqrt(mod$VCV[1,1]), 
                        mod$res_data[1,4],
                        pt(abs(Y_mn/sqrt(mod$VCV[1,1])), df = mod$res_data[1,4], lower.tail = FALSE) * 2)
  
  # Formatting iteration data for outputting
  if(!is.null(plot_dat)){
    if (!is.null(signal_track)) {
      colnames(plot_dat) <- c("Iter", "Conv_check", "FP", "TN", "TP", 
                              "FN", "Total_Disc", "Sum_gamma", "sigma2", "Pred_err", "Conv", "alpha")
      
    } else {
      if (!is.null(report_pred)) {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc", 
                                "Sum_gamma", "sigma2", "Pred_err", "Conv", "alpha")
      } else {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc", 
                                "Sum_gamma", "sigma2", "Conv", "alpha")
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
    cat("Warning: loop completely recycled back to beta=gamma=0 twice. Optimization failed.\n")
  }
  if (try2 == 3) {
    conv <- 3
    cat("Warning: Concentration on 1 variable.\n")
  }
  
  M_step <- LR_update
  Seq_test <- NULL
  
  full_res <- list(beta_ast_hat = beta_ast_hat, beta_hat = beta_hat, beta_hat_var = beta_hat_var, 
                   gamma_hat = gamma, E_step = E_step, Calb_mod = mod, count = count, plot_dat = plot_dat, 
                   Seq_test = Seq_test, M_step = M_step, sigma2_est = mod$sigma2_est, conv = conv, 
                   W_ast = W_ast, W_ast_var = W_ast_var, Y = Y, X = X, 
                   update_order = update_order, beta_start = beta_start, X_mean = X_mean)
  
  # Plotting iteration results if plot_ind=TRUE and either eta_i or test
  # data is given
  if (plot_ind & !is.null(plot_dat)) {
    if (is.null(report_pred)) {
      cat("Warning: cannot plot without eta_i or test data.\n")
    } else {
      plot_probe_func(full_res, test_plot = !is.null(Y_test), alpha = alpha, signal = signal)
    }
  }
  
  return(full_res)
}



