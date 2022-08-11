lr_cpp_func <- function(Y, Z, X = NULL, sigma2) {
  
  N <- length(Y)
  if (!is.null(X)) {
    p <- ncol(X)
    LRcpp <- LM_w_COVS_by_col(Y, Z, as.matrix(X), sigma2)
  } else {
    p <- 0
    LRcpp <- LM_by_col(Y, Z, sigma2)
  }
  t_val <- LRcpp$Coefficients[, 2]/LRcpp$StdErr[, 2]
  
  ret <- list(coef = LRcpp$Coefficients, obs_SE = LRcpp$StdErr, T_val = t_val)
  return(ret)
}