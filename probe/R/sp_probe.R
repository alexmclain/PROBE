#' @title Spatial PaRtitiOned empirical Bayes Ecm (SP-PROBE) algorithm for sparse high-dimensional linear models.
#' @description A wrapper function for the spatially-informed all-at-once variant of the PROBE algorithm.
#' Spatial coordinates are used to borrow strength across neighboring locations when estimating the local
#' false discovery rate (lfdr) in the E-step, increasing power in structured (e.g., imaging) settings.
#' @usage sp_probe(Y, X, coords, Z = NULL, ep = 0.1, maxit = 10000,
#' h_spatial = NULL, Dfull_grid = NULL, h_z = NULL, radius = NULL,
#' null_frac = 0.70, lambda = 0.1, min_neighbors = 10, lambda_tau = 0.1,
#' Y_test = NULL, X_test = NULL, Z_test = NULL, verbose = FALSE,
#' signal = NULL, eta_i = NULL, alpha = 0.05, plot_ind = FALSE, adj = 5)
#' @param Y The outcome variable.
#' @param X An \code{n x M} matrix of sparse predictor variables.
#' @param coords An \code{M x d} matrix of spatial coordinates for each predictor (column of \code{X}),
#' where \code{d} is the number of spatial dimensions (typically 2 for imaging data).
#' @param Z (optional) An \code{n x p} matrix or dataframe of additional predictors not subjected to the sparsity assumption.
#' @param ep Value against which to compare the convergence criterion (default = 0.1).
#' @param maxit Maximum number of iterations the algorithm will run for (default = 10000).
#' @param h_spatial Spatial Gaussian kernel bandwidth controlling the influence of neighboring locations
#' in the local lfdr estimation. If \code{NULL}, estimated as the median pairwise distance divided by 6.
#' @param Dfull_grid (optional) Precomputed full \code{M x M} pairwise distance matrix
#' (e.g., \code{fields::rdist(coords)}). Computed internally from \code{coords} if \code{NULL}.
#' @param h_z KDE bandwidth in the test statistic space. If \code{NULL}, uses a weighted Silverman
#' rule of thumb computed locally within each neighborhood.
#' @param radius Neighbor radius in the same units as \code{coords}. Predictors within this distance
#' of each location are pooled for local estimation. If \code{NULL}, defaults to
#' \code{max(3, min(coord_leng) / 8)}.
#' @param null_frac Fraction of central (near-zero) test statistics used to estimate the local null
#' distribution (default = 0.70).
#' @param lambda Storey tuning parameter for estimating the local proportion of nulls \eqn{\pi_0}
#' (default = 0.1).
#' @param min_neighbors Minimum number of neighbors required for local estimation. If a location has
#' fewer neighbors than this threshold, the full set of predictors is used as a fallback (default = 10).
#' @param lambda_tau Ridge-type penalty applied to the local null variance estimate to prevent
#' instability in sparse neighborhoods (default = 0.1).
#' @param Y_test (optional) Test outcome data, used for plotting purposes only (does not impact results).
#' @param X_test (optional) Test predictor data, used for plotting purposes only (does not impact results).
#' @param Z_test (optional) Test covariate data, used for plotting purposes only (does not impact results).
#' @param verbose Logical; whether to print algorithm iteration progress and summary quantities (default = FALSE).
#' @param signal (optional) A vector of indices of the true non-null coefficients. Used to compute true
#' and false discovery rates by iteration for simulated data. Used for plotting purposes only (does not impact results).
#' @param eta_i (optional) A vector of the true signal values. Used to compute MSE by iteration for
#' simulated data. Used for plotting purposes only (does not impact results).
#' @param alpha Significance level for hypothesis testing (default = 0.05).
#' @param plot_ind Logical; whether to output plots of algorithm results and progress (default = FALSE).
#' @param adj Bandwidth multiplier referenced in the bandwidth adequacy warning message (default = 5).
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
#' @references \itemize{ \item McLain, AC, A Zgodic, H Bondell (2025). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. \emph{Computational Statistics and Data Analysis} 207, 108146.
#' }
#' @examples
#' ### Example
#' data(Sim_data)
#' data(Sim_data_test)
#' 
#' ## Get training and test data
#' Y <- Sim_data$Y
#' X <- Sim_data$X
#' Y_test <- Sim_data_test$Y_test
#' X_test <- Sim_data_test$X_test
#' 
#' ## Get true values
#' signal <- Sim_data$signal
#' beta_tr <- Sim_data$beta_tr
#' sigma2_tr <- Sim_data$sigma2_tr
#' 
#' 
#' # Run the analysis. Y_test and X_test are included for plotting purposes only
#' full_res <- probe( Y = Y, X = X, Y_test = Y_test, 
#' X_test = X_test, alpha = 0.05, plot_ind = TRUE, adj = 10)
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
#' beta_ast_est <- full_res$beta_ast_hat
#' sqrt(mean((beta_ast_est - beta_tr)^2))
#' 
#' # Posterior expectation of gamma by true
#' gamma_est <- full_res$E_step$gamma
#' sum(gamma_est)
#' sum(gamma_est[beta_tr>0])
#' 
#' ### Example with additional covariate data Z (not subjected to the sparsity assumption)
#' data(Sim_data_cov)
#' 
#' # Calculating the true signal (the impact of X only)
#' eta_i <- apply(t(Sim_data_cov$X)*Sim_data_cov$beta_tr,2,sum) 
#  # Run the analysis. eta_i (true signal) and signal are included for plotting purposes only.
#' full_res <- probe( Y = Sim_data_cov$Y, X = Sim_data_cov$X, 
#'                    Z = Sim_data_cov$Z, alpha = 0.05, 
#'                    plot_ind = TRUE, signal = Sim_data_cov$signal, eta_i  = eta_i)
#'                    
#' # Final estimates of the impact of X versus the true values:
#' data.frame(true_values = Sim_data_cov$beta_Z_tr, full_res$Calb_mod$res_data[-2,])
#' 
#' # Compare to a standard linear model of X on Y:
#' summary(lm(Y~Sim_data_cov$Z$Cont_cov + Sim_data_cov$Z$Binary_cov))$coefficients
#' 
#' 
#' @export
sp_probe <- function(Y, X, coords, Z = NULL, ep = 0.1, maxit = 10000, 
                     h_spatial = NULL,    # spatial kernel bandwidth (if NULL use heuristic)
                     Dfull_grid  = NULL,  # Dfull_grid <- fields::rdist(coords) only required when h_spatial = NULL
                     h_z = NULL,          # KDE bandwidth in z (if NULL use Silverman-ish)
                     radius = NULL,       # neighbor radius in grid units (if NULL use 6)
                     null_frac = 0.70,    # fraction of central Z to treat as local null 
                     lambda = 0.1,        # Storey tuning
                     min_neighbors = 10,  # fallback if too few neighbors
                     lambda_tau = 0.1,     # Penalty for estimate of sigma
                  Y_test = NULL, X_test = NULL, Z_test = NULL,
                  verbose = FALSE, signal = NULL, eta_i = NULL, 
                  alpha = 0.05, plot_ind = FALSE) {
  
  if(!is.null(Z)){
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
    }
  }
  
  sp_probe_func(Y = Y, X = X, coords = coords, Z = Z, alpha = alpha, verbose = verbose, 
             h_spatial = h_spatial, h_z = h_z, radius = radius, null_frac = null_frac,
             lambda = lambda, min_neighbors = min_neighbors,
             lambda_tau = lambda_tau, signal = signal, maxit = maxit, eta_i = eta_i,
             ep = ep, plot_ind = plot_ind, Y_test = Y_test, 
             X_test = X_test, Z_test = Z_test)
}


sp_probe_func <- function(Y, X, coords, Z = NULL, alpha, verbose = TRUE, signal, maxit = 1000, 
                       h_spatial = NULL,   
                       h_z = NULL,         
                       radius = NULL,      
                       null_frac = 0.70,   
                       lambda = 0.5,       
                       min_neighbors = 10, 
                       lambda_tau = 0.1,   
                       eta_i = NULL, ep = 0.1, plot_ind = FALSE, 
                       Y_test = NULL, X_test = NULL, Z_test = NULL){
  
  ##### Setting initial values and initializing outputs ####
  M <- dim(X)[2]
  N <- dim(X)[1]
  p <- 0
  if (!is.null(Z)) {
    p <- dim(Z)[2]
    beta_Z <- rep(0, p)
  }
  
  Dfull_grid <- fields::rdist(coords)
  # Lengths of coordinates
  coord_leng <- apply(coords, 2, function(x){length(unique(x))})
  
  # defaults
  # spatial bandwidth heuristic: use half the median neighbor distance among grid-adjacent pixels
  if (is.null(h_spatial)) {
    # median of all pairwise distances (but scaled down)
    h_spatial <- median(Dfull_grid[upper.tri(Dfull_grid)]) / 6
    if (is.na(h_spatial) || h_spatial <= 0) h_spatial <- 3
  }
  
  if (is.null(radius)) radius <- max(3, min(coord_leng) / 8)  # reasonable default
  # precompute neighbor lists (Euclidean radius)
  neighbor_list <- vector("list", N)
  for (i in seq_len(N)) {
    neighbor_list[[i]] <- which(Dfull_grid[,i] <= radius)
  }
  H <- NULL
  if(!is.null(signal)){
    sig_vec <- rep(0,p)
    sig_vec[signal] <- 1
    H <- array(sig_vec, dim = coord_leng)
  }
  
  beta_t <- beta_var <- rep(0,M)
  gamma <- beta_t + 1
  W_ast <- rep(0, N)
  W_ast_var <- W_ast + 1
  sigma2 <- var(Y)
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
      LR_update <- lr_cpp_func(Y, X, Z, sigma2)
      beta_t_new <- c(LR_update$coef[,2])
      # Performing the damping step
      beta_t   <- beta_t*(1-fact) + beta_t_new*fact
      beta_var <- beta_var_old*(1-fact) + (LR_update$obs_SE[,2])^2*fact
    }else {
      LR_update <- m_step_cpp_func(Y, X, Z, W_ast, 
                                   W_ast_var, gamma, 
                                   beta_t, X_2, sigma2) 
      
      beta_t_new <- c(LR_update$coef[,2])
      # Performing the damping step
      beta_t   <- beta_t*(1-fact) + beta_t_new*fact
      beta_var <- 1/((1/beta_var_old)*(1-fact) + (1/LR_update$obs_SE[,2])^2*fact)
    }
    
    
    # Calculating the expectations of the gammas
    E_step <- e_step_func_sp(beta_t, beta_var, df = N - 2, 
                             Dfull_grid = Dfull_grid,
                             neighbor_list = neighbor_list,
                             h_spatial = h_spatial, h_z = h_z, 
                             radius = radius, 
                             null_frac = null_frac,
                             lambda = lambda, 
                             min_neighbors = min_neighbors,
                             lambda_tau = lambda_tau,  
                             coord_leng = coord_leng,
                             H = H) 
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
      mod <- m_step_regression(Y, W_ast, W_ast_var + W_ast^2, Z, a = -3/2, Int = TRUE) 
      Y_pred <- mod$Y_pred
      sigma2 <- mod$sigma2_est 
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
            message("Iteration=", count, " Number of discoveries (using lfdr)=", 
                    disc, " Sum(gamma)=", round(sum(E_step$gamma), 1), " MSPE(test)=", 
                    report_pred, " Convergence Crit=", CC_round, "\n")
          } else {
            message("Iteration=", count, " Number of discoveries (using lfdr)=", 
                    disc, " Sum(gamma)=", round(sum(E_step$gamma), 1), " Convergence Crit=", 
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
            too narrow, consider raising 'adj'. Current adj=", h_z))
  }
  
  
  full_res <- list(beta_ast_hat = beta_ast_hat, beta_hat = beta_hat, 
                   beta_hat_var = beta_hat_var, 
                   gamma_hat = gamma, E_step = E_step, 
                   Calb_mod = mod, count = count, plot_dat = plot_dat, 
                   M_step = LR_update, sigma2_est = mod$sigma2_est, conv = conv, 
                   W_ast = W_ast, W_ast_var = W_ast_var, Y = Y, X = X, Z = Z)
  
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



