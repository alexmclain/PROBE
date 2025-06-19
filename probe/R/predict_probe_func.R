#' @title Obtaining predictions, confidence intervals and prediction intervals from probe
#' @description A function providing predictions, along with \eqn{(1-\alpha)*100\%} credible, and prediction intervals for new observations.
#'
#'
#' @usage predict_probe_func(res, X, Z = NULL, alpha = 0.05, X_2 = NULL)
#' @param res The results from the probe function.
#' @param X A matrix containing the predictors on which to apply the probe algorithm
#' @param Z (optional) A matrix or dataframe of predictors not subjected to the sparsity assumption to account for.
#' @param alpha significance level for (\eqn{100(1-\alpha)\%}) credible and prediction intervals.
#' @param X_2 (optional) Square of X matrix.
#' @return A dataframe with predictions, credible intervals, and prediction intervals for each new observation.
#' @references \itemize{ \item McLain, AC, A Zgodic, H Bondell (2025). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. \textit{Computational Statistics and Data Analysis} 207, 108146.
#' \item Zgodic, A., Bai, R., Zhang, J., Wang, Y., Rorden, C., & McLain, A. (2023). Quantifying predictive uncertainty of aphasia severity in stroke patients with sparse heteroscedastic Bayesian high-dimensional regression. arXiv preprint arXiv:2309.08783.}
#' @examples
#' ### Example
#' #not run
#' # pred_res <- predict_probe_func(full_res, X = X_test, Z = NULL, alpha = alpha)
#' # head(pred_res)
#' 
#' @export
predict_probe_func <- function(res, X, Z = NULL, alpha = 0.05, 
                                 X_2 = NULL) {
  
  if(is.data.frame(X)){stop("Error: X must be a matrix or vector.")}
  M <- length(res$beta_hat)
  E_step <- res$E_step
  mod <- res$Calb_mod
  sigma2_est <- mod$sigma2_est
  coef_est <- mod$coef
  VCV <- mod$VCV
  if (!is.matrix(X)){
    if(!is.null(dim(X))[1]){
      N_test <- prod(dim(X))/M
    }else{
      N_test <- 1
    }
    X <- matrix(X,N_test,M)
    X_2 <- NULL
  }
  if(!is.null(res$X_mean)){
    X <- t(t(X) - res$X_mean)
  }
  
  if (is.null(X_2)) X_2 <- X * X
  X <- as.matrix(X)
  X_2 <- as.matrix(X_2)
  
  W_W2_update <- m_update_func(X, X_2, E_step$beta_tilde, E_step$gamma, E_step$beta_tilde_var)
  W_ast <- W_W2_update$W_ast
  W_ast_var <- W_W2_update$W_ast_var
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
  }
  mod_mat <- cbind(1, W_ast, Z)
  pred_mean <- mod_mat %*% coef_est
  Var_train <- diag(mod_mat %*% VCV %*% t(mod_mat)) + 
    (VCV[2, 2] + coef_est[2]^2) * W_ast_var
  
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
