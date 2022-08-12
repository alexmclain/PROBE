#' @title Obtaining predictions, confidence intervals and prediction intervals from probe
#' @description A function providing predictions, along with \eqn{(1-\alpha)*100\%} credible, and prediction intervals for new observations.
#'
#'
#' @usage predict_probe_func(res, Z, X, alpha, Z_2)
#' @param res The results from the probe function.
#' @param Z A matrix containing the predictors on which to apply the probe algorithm
#' @param X (optional) A matrix or dataframe of predictors not subjected to the sparsity assumption to account for.
#' @param alpha Type I error; significance level
#' @param Z_2 (optional) Square of Z matrix.
#' @return A dataframe with predictions, confidence intervals, and prediction intervals for each new observation.
#' @examples
#' ### Example
#' data(Sim_data)
#' attach(Sim_data)
#' alpha <- 0.05
#' plot_ind <- TRUE
#' 
#' # Run the analysis. Y_test and Z_test are included for plotting purposes only
#' full_res <- probe( Y = Y, Z = Z, alpha = alpha, plot_ind = plot_ind, Y_test = Y_test, Z_test = Z_test)
#' 
#' # Predicting for test data
#' pred_res <- predict_probe_func(full_res, Z = Z_test, X = NULL, alpha = alpha)
#' head(pred_res)
#' 
predict_probe_func <- function(res, Z, X = NULL, alpha = 0.05, 
                                 Z_2 = NULL) {
  E_step <- res$E_step
  mod <- res$Calb_mod
  sigma2_est <- res$Calb_mod$sigma2_est
  if (is.null(Z_2)) Z_2 <- Z * Z
  
  W_W2_update <- m_update_func(Z, Z_2, E_step$beta_tilde, E_step$delta, E_step$beta_tilde_var)
  W_ast <- W_W2_update$W_ast
  W_ast_var <- W_W2_update$W_ast_var
  if (!is.null(X)) {
    X <- as.matrix(X)
  }
  mod_mat <- cbind(1, W_ast, X)
  pred_mean <- mod_mat %*% mod$coef
  Var_train <- diag(mod_mat %*% mod$VCV %*% t(mod_mat)) + 
    (mod$VCV[2, 2] + mod$coef[2]^2) * W_ast_var
  
  pred_res <- data.frame(Pred = pred_mean)
  CI_train <- PI_train <- NULL
  if (!is.null(alpha)) {
    pred_res$CI_L <- pred_mean - qnorm(1 - alpha/2) * sqrt(Var_train)
    pred_res$CI_U <- pred_mean + qnorm(1 - alpha/2) * sqrt(Var_train)
    pred_res$PI_L <- pred_mean - qnorm(1 - alpha/2) * sqrt(Var_train + sigma2_est)
    pred_res$PI_U <- pred_mean + qnorm(1 - alpha/2) * sqrt(Var_train + sigma2_est)
  }
  pred_res$Var <- Var_train
  pred_res$Var_pred <- Var_train + sigma2_est
  
  return(pred_res)
}
