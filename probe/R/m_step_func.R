#' @title Function for fitting the initial part of the M-step
#' @description A wrapper function providing the quantities related to the M-step for \eqn{\alpha_0} and \eqn{\sigma^2}.
#' @param Y A matrix containing the outcome \code{Y}
#' @param W Quantity \eqn{E(W_0)} as outlined in citation, output from W_update_fun
#' @param W2 Quantity \eqn{E(W^2_0)} as outlined in citation, output from W_update_fun
#' @param Z A matrix or dataframe of other predictors to account for
#' @param a (optional) parameter for changing the hyperparameter \eqn{a} (default, \eqn{a=-3/2} uses \eqn{n-2} as denominator for MAP of \eqn{\sigma^2})
#' @param Int (optional) Logical - should an intercept be used?
#' @return A list including 
#' 
#' \code{coef} the MAP estimates of the \eqn{\alpha_0} parameters
#' \code{sigma2_est} the MAP estimate of \eqn{\sigma^2}
#' \code{VCV} posterior variance covariance matrix of \eqn{\alpha_0}, 
#' \code{res_data} dataframe containing MAP estimates, posterior variances, t-test statistics and associated p-values for \eqn{\alpha_0}
#' @references 
#' McLain, A. C., Zgodic, A., & Bondell, H. (2022). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2209.08139.
#' @examples
#' #not run
#' #mod <- m_step_regression(Y, W_ast, W_ast_var + W_ast^2, Z)
m_step_regression <- function(Y, W, W2, Z = NULL, a = -3/2, Int = TRUE) {
  
  N <- length(Y)
  if(!is.null(W)){
    if(Int){ Wmat <- cbind(1, W, Z) }
    if(!Int & !is.null(Z) ){ Wmat <- cbind(W, Z) }
    if(!Int & is.null(Z) ){ Wmat <- matrix(W) }
    W_col <- c(1 + 1*I(Int))
    
    WWpr  <- t(Wmat) %*% Wmat
    WWpr[W_col, W_col] <- sum(W2)
    
    df       <- N - ncol(Wmat) + 2*a - 3
    
    if (det(WWpr) != 0) {
      WWpr_inv <- solve(WWpr)
      beta_w <- WWpr_inv %*% t(Wmat) %*% Y
      Y_pred <- Wmat %*% beta_w
      hat <- diag(Wmat %*% WWpr_inv %*% t(Wmat))
      resid <- (Y - Y_pred)
      RSS <- sum(resid^2) + beta_w[W_col]^2*(sum(W2) - sum(W^2))
      VCV <- RSS/(df) * t(WWpr_inv) %*% (t(Wmat) %*% Wmat) %*% WWpr_inv
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
  list(coef = beta_w, sigma2_est = RSS/df, RSS = RSS, 
       Y_pred = Y_pred, hat = hat, resid = resid, 
       VCV = VCV, res_data = res_data)
  
}

m_step_cpp_func <- function(Y, X, Z = NULL, W_ast, W_ast_var, gamma, 
                            beta_vec, X_2, sigma2) {
  
  N <- length(Y)
  if (!is.null(Z)) {
    LRcpp <- PROBE_cpp0_5_6_covs(Y, X, W_ast, W_ast_var, gamma, 
                                 beta_vec, X_2, sigma2, as.matrix(Z))
  } else {
    
    LRcpp <- PROBE_cpp0_5_6(Y, X, W_ast, W_ast_var, gamma, 
                            beta_vec, X_2, sigma2)
  }
  
  ret <- list(coef = LRcpp$Coefficients, obs_SE = LRcpp$StdErr)
  
  return(ret)
}

m_step_cpp_one <- function(Y, X, W_ast, W_ast_var, gamma, 
                             beta_t, X_2, sigma2, update_order) {
  
  N <- length(Y)
  LRcpp <- PROBE_one_cpp(Y, X, W_ast, W_ast_var, gamma, 
                           beta_t, X_2, sigma2, update_order)
  
  LRcpp
}


lr_cpp_func <- function(Y, X, Z = NULL, sigma2) {
  
  N <- length(Y)
  if (!is.null(Z)) {
    LRcpp <- LM_w_COVS_by_col(Y, X, as.matrix(Z), sigma2)
  } else {
    LRcpp <- LM_by_col(Y, X, sigma2)
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



