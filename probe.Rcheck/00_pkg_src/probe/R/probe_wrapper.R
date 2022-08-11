#' @title Fitting PaRtitiOned empirical Bayes Ecm (PROBE) algorithm to sparse high-dimensional linear models.
#' @description A wrapper function for the main PROBE algorithm function.  The R package is a work in progress.
#' @usage probe(Y, Z, X, alpha, ep, maxit, Y_test, Z_test, X_test, verbose, signal, eta_i, plot_ind)
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
#' * \code{beta_ast_hat} MAP estimates of the regression coefficients (\eqn{\beta^\ast}),
#' * \code{beta_hat, beta_hat_var} MAP estimates of the posterior expectation (beta_hat) and variance (beta_hat_var) of the prior mean (\eqn{\beta}) of the regression coefficients assuming \eqn{\gamma=1}, 
#' * \code{gamma_hat} the posterior expectation of the latent \eqn{\gamma} variables, 
#' * \code{sigma2_est} MAP estimate of the residual variance, 
#' * \code{E_step} full results of the final E_step, 
#' * \code{Calb_mod} results of first (\eqn{\alpha_0}) part of the M-step,  
#' * \code{count} the total number of iterations before convergence. 
#' 
#' @seealso predict_probe_func to obtain predictions, confidence intervals and prediction intervals from PROBE.
#' @examples
#' ### Example
#' data(Sim_data)
#' attach(Sim_data)
#' alpha <- 0.05
#' plot_ind <- TRUE
#' 
#' # Run the analysis. Y_test and Z_test are included for plotting purposes only
#' full_res <- probe( Y = Y, Z = Z, alpha = alpha, Y_test = Y_test, Z_test = Z_test, plot_ind = plot_ind)
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
#' full_res <- probe( Y = Y, Z = Z, X = X, alpha = alpha, signal = signal, eta_i  = eta_i, plot_ind = plot_ind)
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
