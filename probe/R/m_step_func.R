#' A wrapper function providing the quantities related to the M-step for \eqn{\alpha_0} and \eqn{\sigma^2}.
#'
#' @param Y A matrix containing the outcome \code{Y}
#' @param W Quantity \eqn{E(W_0)} as outlined in citation, output from W_update_fun
#' @param W2 Quantity \eqn{E(W^2_0)} as outlined in citation, output from W_update_fun
#' @param X A matrix or dataframe of other predictors to account for
#' @param a (optional) parameter for changing the hyperparameter \eqn{a} (default, \eqn{a=-3/2} uses \eqn{n-2} as denominator for MAP of \eqn{\sigma^2})
#' @param Int (optional) Logical - should an intercept be used?
#' @return A list including 
#' 
#' \code{coef} the MAP estimates of the \eqn{\alpha_0} parameters
#' \code{sigma2_est} the MAP estimate of \eqn{\sigma^2}
#' \code{VCV} posterior variance covariance matrix of \eqn{\alpha_0}, 
#' \code{res_data} dataframe containing MAP estimates, posterior variances, t-test statistics and associated p-values for \eqn{\alpha_0}
#' @examples
#' #not run
#' #mod <- m_step_regression(Y, W_ast, W_ast_var + W_ast^2, X)
m_step_regression <- function(Y, W, W2, X, a = -3/2, Int = TRUE) {
  
  N <- length(Y)
  if(!is.null(W)){
    if(Int){ Wmat <- cbind(1, W, X) }
    if(!Int & !is.null(X) ){ Wmat <- cbind(W, X) }
    if(!Int & is.null(X) ){ Wmat <- matrix(W) }
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
      RSS <- sum(resid^2)
      SST <- sum((Y - mean(Y))^2)
      R_squared <- 1 - RSS/SST
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
      SST <- sum((Y - mean(Y))^2)
      R_squared <- 1 - RSS/SST
      VCV <- vcov(final_mod)
    }
    
    Std_Err  <- sqrt(diag(VCV))
    T_vals   <- beta_w/Std_Err
  }
  if(is.null(W)){
    SST <- sum((Y - mean(Y))^2)
    if(is.null(X)){
      fail_lm <- lm(Y~1)
      coef_fail <- summary(fail_lm)$coefficients
      Y_pred <- as.numeric(fail_lm$fitted.values)
      hat <- as.numeric(influence(fail_lm)$hat)
      resid <- as.numeric(fail_lm$residuals)
      RSS <- sum(resid^2)
      R_squared <- 1 - RSS/SST
      df <- N- 1 - + 2*a - 3
      beta_w <- c(coef_fail[1], 0)
      VCV <- matrix(0,2,2)
      VCV[1,1] <- vcov(fail_lm)
      Std_Err  <- c(sqrt(vcov(fail_lm)), 0)
      T_vals   <- c(beta_w[1]/Std_Err[1], 0)
    }else{
      Wmat <- cbind(1, X)
      fail_lm <- lm(Y~Wmat)
      coef_fail <- summary(fail_lm)$coefficients
      Y_pred <- as.numeric(fail_lm$fitted.values)
      hat <- as.numeric(influence(fail_lm)$hat)
      resid <- as.numeric(fail_lm$residuals)
      RSS <- sum(resid^2)
      R_squared <- 1 - RSS/SST
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
  if(Int){row.names(res_data) <- c("Intercept", "W", colnames(X))}
  if(!Int){row.names(res_data) <- c("W", colnames(X))}
  list(coef = beta_w, sigma2_est = RSS/df, RSS = RSS, R_squared = 
         R_squared, Y_pred = Y_pred, hat = hat, resid = resid, 
       VCV = VCV, res_data = res_data)
  
}