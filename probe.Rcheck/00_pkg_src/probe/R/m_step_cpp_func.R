m_step_cpp_func <- function(Y, Z, X, W_ast, W_ast_var, delta, 
                            beta_vec, Z_2, sigma2) {

  N <- length(Y)
  if (!is.null(X)) {
    p <- ncol(X)
    LRcpp <- PROBE_cpp0_5_6_covs(Y, Z, W_ast, W_ast_var, delta, 
                                   beta_vec, Z_2, sigma2, as.matrix(X))
  } else {
    p <- 0
    LRcpp <- PROBE_cpp0_5_6(Y, Z, W_ast, W_ast_var, delta, 
                              beta_vec, Z_2, sigma2)
  }
  t_val <- LRcpp$T_statistics
  
  ret <- list(coef = LRcpp$Coefficients, obs_SE = LRcpp$StdErr, T_val = t_val)
  
  return(ret)
}