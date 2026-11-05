#' @title Spatial PaRtitiOned empirical Bayes Ecm (SP-PROBE) algorithm for sparse high-dimensional linear models.
#' @description A wrapper function for the spatially-informed all-at-once variant of the PROBE algorithm.
#' Spatial coordinates are used to borrow strength across neighboring locations when estimating the local
#' false discovery rate (lfdr) in the E-step, increasing power in structured (e.g., imaging) settings.
#' @usage sp_probe(Y, X, coords, Z = NULL, maxit = 10000, L = 50,
#' h_spatial = NULL, Dfull_grid = NULL, h_z = 1, radius = NULL,
#' null_frac = 0.70, lambda = 0.1, min_neighbors = 10, lambda_tau = 0.001,
#' Y_test = NULL, X_test = NULL, Z_test = NULL, verbose = FALSE,
#' signal = NULL, eta_i = NULL, alpha = 0.05, plot_ind = FALSE,
#' method = c("kde", "cnf"), estimate_sigma = TRUE, sigma_null = 1,
#' cnf_K = 10, cnf_hidden_dim = 64, cnf_n_hidden = 2,
#' cnf_n_epochs = 1000, cnf_lr = 1e-3, cnf_batch_size = 256,
#' monotone = FALSE, monotone_step = 0.1)
#' @param Y The outcome variable.
#' @param X An \code{n x M} matrix of sparse predictor variables.
#' @param coords An \code{M x d} matrix of spatial coordinates for each predictor (column of \code{X}),
#' where \code{d} is the number of spatial dimensions (typically 2 for imaging data).
#' @param Z (optional) An \code{n x p} matrix or dataframe of additional predictors not subjected to the sparsity assumption.
#' @param maxit Maximum number of iterations the algorithm will run for (default = 10000).
#' @param L Patience window: number of consecutive iterations without improvement in \eqn{\hat{\sigma}^2}
#' before the algorithm stops and returns the coefficients from the iteration with the minimum
#' \eqn{\hat{\sigma}^2} (default = 50). Only used when \code{conv_method = "sigma2"}.
#' @param conv_method Convergence criterion to use. \code{"sigma2"} (default) stops after \code{L}
#' consecutive iterations without improvement in \eqn{\hat{\sigma}^2} and returns the state from
#' the iteration with the minimum \eqn{\hat{\sigma}^2}. \code{"chi2"} stops when a \eqn{\chi^2}
#' statistic measuring the change in \eqn{\tilde{W}} between consecutive iterations falls below
#' \code{ep} for two successive iterations, and returns the state at convergence.
#' @param ep Convergence threshold for \code{conv_method = "chi2"} (default = 0.1). The
#' \eqn{\chi^2} statistic is \code{pchisq(max((W_old - W)^2 / Var_W) / log(N), df = 1)};
#' convergence is declared when this quantity and its value from the previous iteration are both
#' below \code{ep}.
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
#' @param monotone Logical; only used when \code{method = "kde"}. If \code{TRUE}, each location's
#' lfdr is made non-increasing in \eqn{|T|} at every E-step, as in PROBE's \code{monotone = TRUE}
#' (see \code{e_step_func_sp}). Considerably slower (default = \code{FALSE}).
#' @param monotone_step Grid spacing in \eqn{t} used when \code{monotone = TRUE} (default = 0.1).
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
#' When \code{TRUE} (and \code{method = "kde"} on a 2D grid), also draws the 3 x 3 spatial E-step
#' diagnostic panel (lfdr, f, f0, pi0, bandwidth, local null SD) for the returned state.
#' @param adj Bandwidth multiplier referenced in the bandwidth adequacy warning message (default = 5).
#' @param method Density estimation method for the marginal distribution in the E-step.
#' \code{"kde"} (default) uses spatially-weighted kernel density estimation via
#' \code{e_step_func_sp}. \code{"cnf"} uses a conditional normalizing flow
#' (neural-network mixture-of-logistics); requires the \code{torch} package and the
#' helper functions \code{train_cnf}, \code{estimate_f_cnf}, and \code{lfdr_threshold}.
#' When \code{method = "cnf"} the CNF model is trained once on the first iteration and
#' reused in all subsequent iterations.
#' @param estimate_sigma Logical; only used when \code{method = "cnf"}. If \code{TRUE}
#' (default), estimates the local null variance via a weighted truncated-normal MLE. If
#' \code{FALSE}, fixes null SD to \code{sigma_null}.
#' @param sigma_null Fixed null SD used when \code{method = "cnf"} and
#' \code{estimate_sigma = FALSE} (default = 1).
#' @param cnf_K Number of mixture-of-logistics components for the CNF (default = 10).
#' @param cnf_hidden_dim Conditioner network width for the CNF (default = 64).
#' @param cnf_n_hidden Number of hidden layers in the CNF conditioner (default = 2).
#' @param cnf_n_epochs Adam training epochs for the CNF (default = 1000).
#' @param cnf_lr Adam learning rate for the CNF (default = 1e-3).
#' @param cnf_batch_size Mini-batch size for CNF training (default = 256).
#' @param cnf_update_rep How often the global CNF model is updated (default = 1, updated fresh at each iteration).
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
#' @importFrom fields rdist
#' @export
sp_probe <- function(Y, X, coords, Z = NULL, maxit = 10000, L = 50,
                     conv_method = c("sigma2", "chi2"),
                     ep = 0.1,
                     h_spatial = NULL,    # spatial kernel bandwidth (if NULL use heuristic)
                     Dfull_grid  = NULL,  # precomputed M x M distance matrix; computed from coords if NULL
                     h_z = 1,             # KDE bandwidth multiplier in z space (method = "kde" only)
                     radius = NULL,       # neighbor radius in coord units (if NULL use heuristic)
                     null_frac = 0.70,    # fraction of central Z to treat as local null
                     lambda = 0.1,        # Storey tuning
                     min_neighbors = 10,  # fallback if too few neighbors
                     lambda_tau = 0.001,  # ridge penalty for null sigma estimate
                     Y_test = NULL, X_test = NULL, Z_test = NULL,
                     verbose = FALSE, signal = NULL, eta_i = NULL,
                     alpha = 0.05, plot_ind = FALSE,
                     method = c("kde", "cnf"),
                     estimate_sigma = TRUE,
                     sigma_null = 1,
                     cnf_K = 10,
                     cnf_hidden_dim = 64,
                     cnf_n_hidden = 2,
                     cnf_n_epochs = 1000,
                     cnf_lr = 1e-3,
                     cnf_batch_size = 256, 
                     cnf_update_rep = 1,
                     monotone = FALSE,
                     monotone_step = 0.1) {

  conv_method <- match.arg(conv_method)
  method      <- match.arg(method)
  if(!is.null(Z)){
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
    }
  }

  sp_probe_func(Y = Y, X = X, coords = coords, Dfull_grid = Dfull_grid, Z = Z, alpha = alpha, verbose = verbose,
                h_spatial = h_spatial, h_z = h_z, radius = radius, null_frac = null_frac,
                lambda = lambda, min_neighbors = min_neighbors,
                lambda_tau = lambda_tau, signal = signal, maxit = maxit, eta_i = eta_i,
                L = L, conv_method = conv_method, ep = ep,
                plot_ind = plot_ind, Y_test = Y_test,
                X_test = X_test, Z_test = Z_test,
                method = method, estimate_sigma = estimate_sigma, sigma_null = sigma_null,
                cnf_K = cnf_K, cnf_hidden_dim = cnf_hidden_dim, cnf_n_hidden = cnf_n_hidden,
                cnf_n_epochs = cnf_n_epochs, cnf_lr = cnf_lr, cnf_batch_size = cnf_batch_size,
                monotone = monotone, monotone_step = monotone_step)
}


sp_probe_func <- function(Y, X, coords,
                          Dfull_grid = NULL,
                          Z = NULL, alpha,
                          verbose = TRUE, signal, maxit = 1000,
                          h_spatial = NULL,
                          h_z = 1,
                          radius = NULL,
                          null_frac = 0.70,
                          lambda = 0.1,
                          min_neighbors = 10,
                          lambda_tau = 0.001,
                          eta_i = NULL, L = 50,
                          conv_method = c("sigma2", "chi2"), ep = 0.1,
                          plot_ind = FALSE,
                          Y_test = NULL, X_test = NULL, Z_test = NULL,
                          method = c("kde", "cnf"),
                          estimate_sigma = TRUE,
                          sigma_null = 1,
                          cnf_K = 10,
                          cnf_hidden_dim = 64,
                          cnf_n_hidden = 2,
                          cnf_n_epochs = 1000,
                          cnf_lr = 1e-3,
                          cnf_batch_size = 256, 
                          cnf_update_rep = 1,
                          monotone = FALSE,
                          monotone_step = 0.1){
  conv_method <- match.arg(conv_method)
  method      <- match.arg(method)
  
  ##### Setting initial values and initializing outputs ####
  M <- dim(X)[2]
  N <- dim(X)[1]
  p <- 0
  if (!is.null(Z)) {
    p <- dim(Z)[2]
    beta_Z <- rep(0, p)
  }
  
  if(verbose){message("Setting up spatial coordinates.\n")}
  if(is.null(Dfull_grid)){Dfull_grid <- fields::rdist(coords)}
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
  neighbor_list <- vector("list", M)
  for (i in seq_len(M)) {
    neighbor_list[[i]] <- which(Dfull_grid[,i] <= radius)
  }
  H <- NULL
  if(!is.null(signal)){
    sig_vec <- rep(0,M)
    sig_vec[signal] <- 1
    H <- array(sig_vec, dim = coord_leng)
  }
  # spatial KDE weights are fixed across iterations, so compute them once
  sp_weights <- NULL
  if (method == "kde") sp_weights <- build_sp_weights(neighbor_list, Dfull_grid, h_spatial, min_neighbors)
  if (verbose) {
    message("spatial_nonparam_kde: radius = ", round(radius, 3),
            ", h_spatial = ", round(h_spatial,3),
            ", h_z = ", round(h_z,2))
    message("On to the algorithm")}
  beta_t <- beta_var <- rep(0,M)
  gamma <- beta_t + 1
  W_ast <- rep(0, N)
  W_ast_var <- W_ast + 1
  sigma2 <- var(Y)
  count <- rcy_ct <- conv_check <- try2 <- 0
  plot_dat <- NULL
  X_2 <- X * X
  sigma2_min <- Inf
  no_improve <- 0L
  best_state <- NULL
  RSS <- NA
  Xt_conv1 <- prev_Xt_conv1 <- 1
  T_vals <- signal_track <- report_pred <- plot_dat <- NULL
  cnf_model <- NULL  # CNF model trained on first iteration and reused thereafter
  
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
    t_verbose <- FALSE
    if ((count %% 25) == 0 & verbose){
      t_verbose = TRUE
    } 
    
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
                             H = H,
                             verbose = t_verbose,
                             sp_weights = sp_weights,
                             monotone = monotone,
                             monotone_step = monotone_step,
                             method = method,
                             coords = coords,
                             estimate_sigma = estimate_sigma,
                             sigma_null = sigma_null,
                             cnf_model = cnf_model,
                             cnf_K = cnf_K, cnf_hidden_dim = cnf_hidden_dim,
                             cnf_n_hidden = cnf_n_hidden, cnf_n_epochs = cnf_n_epochs,
                             cnf_lr = cnf_lr, cnf_batch_size = cnf_batch_size)
    # Cache the fitted CNF model so it is reused in subsequent iterations
    # The model is only reused if count %% cnf_update_rep != 0, so 
    # if cnf_update_rep = 1 the model is never reused
    # if cnf_update_rep = 5 the model is re-updated every 5 iterations
    
    if (method == "cnf" && !is.null(E_step$cnf_model)) {
      if((count %% cnf_update_rep) == 0 | count <= 10){
        cnf_model <- NULL
      }else{cnf_model <- E_step$cnf_model}
    }
    t_verbose <- FALSE
    gamma <- E_step$gamma
    
    if (sum(gamma) == 0) {
      if(try2 == 0){
        warning("algorithm recycled back to null gamma values.\n 
            Trying again with different starting values. \n")
        rcy_ct <- count
        beta_t <- rep(1e-5,M)
        gamma <- rep(1,M)
        T_vals <- NULL
        count <- try2 <- 1
        no_improve <- 0L; sigma2_min <- Inf; best_state <- NULL
        cnf_model <- NULL
        Xt_conv1 <- prev_Xt_conv1 <- 0
      } else {
        E_step$beta_tilde <-  beta_t
        E_step$gamma <-  gamma_old
        E_step$beta_tilde_var <- beta_var
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
      RSS    <- mod$RSS
    }
    
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
    
    # Always track the state from the iteration with minimum sigma2.
    # Used as the return value for conv_method="sigma2" and as a fallback
    # for conv_method="chi2" when the algorithm does not converge.
    if (var(W_ast) > 0) {
      if (sigma2 < sigma2_min) {
        sigma2_min <- sigma2
        no_improve <- 0L
        best_state <- list(beta_t       = beta_t,    beta_var     = beta_var,
                           gamma        = gamma,      mod          = mod,
                           W_ast        = W_ast,      W_ast_var    = W_ast_var,
                           E_step       = E_step,     LR_update    = LR_update,
                           count        = count,
                           report_pred  = report_pred,
                           MTR_res      = MTR_res,
                           signal_track = signal_track)
      } else {
        no_improve <- no_improve + 1L
        if (conv_method == "sigma2" && no_improve >= L) conv_check <- conv_check + 1L
      }
    }
    
    # Chi-square convergence check (conv_method = "chi2" only).
    # Declares convergence when the chi2 statistic on the change in W_ast
    # is below ep for two consecutive iterations.
    if (conv_method == "chi2" && count > 1) {
      if (length(c(W_ast_old[W_ast2_old > 0], W_ast[W_ast2_old > 0])) > 0) {
        Xt_conv1 <- pchisq(max((W_ast_old[W_ast2_old > 0] - W_ast[W_ast2_old > 0])^2 /
                                 W_ast2_old[W_ast2_old > 0]) / log(N), df = 1)
      }
      if (Xt_conv1 < ep) { conv_check <- conv_check + 1L }
      prev_Xt_conv1 <- Xt_conv1
      no_improve <- Xt_conv1
    }
    
    # Outputting results if verbose=TRUE
    if (verbose & ((count %% 25) == 0 | conv_check == 1)) {
      if (conv_method == "sigma2") {
        report_pred_sp(no_improve, L, best_state$signal_track, best_state$MTR_res,
                       best_state$report_pred, best_state$count, best_state$E_step)
      } else {
        report_pred_sp(round(Xt_conv1, 4), ep, signal_track, MTR_res,
                       report_pred, count, E_step, label = "chi2_stat")
      }
    }
    # Storing iteration data
    plot_dat <- rbind(plot_dat, c(rcy_ct+count, conv_check, signal_track,
                                  sum(MTR_res$BH_res$LFDR), sum(gamma),
                                  sigma2, report_pred, no_improve, RSS))
  }
  
  # Restore state.
  # conv_method="sigma2": always use the min-sigma2 best_state.
  # conv_method="chi2", converged: use the current iteration's state.
  # conv_method="chi2", not converged (hit maxit): fall back to best_state.
  if (conv_method == "sigma2" || conv_check == 0) {
    if (is.null(best_state)) stop("sp_probe: sigma2 was never estimated; check that var(W_ast) > 0.")
    beta_hat     <- best_state$beta_t
    beta_hat_var <- best_state$beta_var
    gamma        <- best_state$gamma
    mod          <- best_state$mod
    W_ast        <- best_state$W_ast
    W_ast_var    <- best_state$W_ast_var
    E_step       <- best_state$E_step
    LR_update    <- best_state$LR_update
    count        <- best_state$count
  } else {
    # chi2 method, converged: current iteration state is the converged state
    beta_hat     <- beta_t
    beta_hat_var <- beta_var
    # gamma, mod, W_ast, W_ast_var, E_step, LR_update, count are already current
  }
  beta_ast_hat <- beta_hat * gamma * mod$coef[2]
  
  
  # Formatting iteration data for outputting
  if (!is.null(plot_dat)) {
    conv_col <- if (conv_method == "sigma2") "Patience" else "chi2_stat"
    if (!is.null(signal_track)) {
      colnames(plot_dat) <- c("Iter", "Conv_check", "FP", "TN", "TP",
                              "FN", "Total_Disc", "Sum_gamma", "sigma2", "Pred_err", conv_col, "RSS")
    } else {
      if (!is.null(report_pred)) {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc",
                                "Sum_gamma", "sigma2", "Pred_err", conv_col, "RSS")
      } else {
        colnames(plot_dat) <- c("Iter", "Conv_check", "Total_Disc",
                                "Sum_gamma", "sigma2", conv_col, "RSS")
      }
    }
    plot_dat <- data.frame(plot_dat)
  }
  
  conv <- 0
  if (conv_check == 0 & try2 < 2) {
    conv <- 1
    warning("convergence criteria not met. Set different convergence criteria or raise maximum number of iterations.")
  }
  if (try2 == 2) {
    conv <- 2
    warning("algorithm recycled back to null values twice. Optimization failed.")
  }
  
  if(E_step$adj_warning != 0){
    message("Results indicate bandwidth in density estimation may be too narrow.\n Consider raising h_z, lambda_tau, or using monotone = TRUE.")
  }
  
  
  full_res <- list(beta_ast_hat = beta_ast_hat, beta_hat = beta_hat, 
                   beta_hat_var = beta_hat_var, 
                   gamma_hat = gamma, E_step = E_step, 
                   Calb_mod = mod, count = count, plot_dat = plot_dat, 
                   M_step = LR_update, sigma2_est = mod$sigma2_est, conv = conv, 
                   W_ast = W_ast, W_ast_var = W_ast_var, Y = Y, X = X, Z = Z)
  
  # Plotting the spatial E-step diagnostics for the returned state
  if (plot_ind && method == "kde" && length(coord_leng) == 2 && !is.null(E_step$EB_res)) {
    plot_e_step_sp(E_step$EB_res, E_step$T_vals, E_step$p_vals, Dfull_grid, coord_leng, H = H)
  }

  # Plotting iteration results if plot_ind=TRUE and either eta_i or test
  # data is given
  if (plot_ind & !is.null(plot_dat)) {
    if (is.null(report_pred)) {
      warning("cannot plot without true mean or test data.\n")
    } else {
      # conv_method="sigma2": drop the L patience iterations run after the
      # returned (min-sigma2) iteration; plot_dat in full_res keeps them.
      plot_res <- full_res
      if (conv_method == "sigma2") {
        plot_res$plot_dat <- plot_dat[plot_dat$Iter <= rcy_ct + count, , drop = FALSE]
      }
      plot_probe_func(plot_res, test_plot = !is.null(Y_test), alpha = alpha, signal = signal)
    }
  }
  
  
  return(full_res)
}



