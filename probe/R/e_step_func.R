e_step_func <- function(beta_t, beta_var, df, adj, lambda = 0.1, monotone=FALSE){
  
  ### Get and clean up T statistics
  T_vals <- beta_t/sqrt(beta_var)
  T_vals[beta_var<=0] <- 0
  T_vals[is.na(T_vals)] <- 0
  T_vals[is.nan(T_vals)] <- 0   
  M <- length(T_vals) 
  p_vals <- pt(abs(T_vals),df=df,lower.tail = FALSE)*2

  ### Estimate the probability of a null
  pi_hat  <- (pi0_func(p_vals,lambda = lambda)$pi0*M)/(M+1)
  ### Estimate the lfdr
  delta <- 1-lfdr_t_func(T = T_vals, pi0 = pi_hat, trunc=TRUE, 
                         monotone = monotone, adj = adj, 
                         df_val = df)$lfdr
  
  ret <- list(beta_tilde = beta_t, beta_tilde_var = beta_var, delta = delta, 
              lfdr = 1 - delta, pi0 = pi_hat, p_vals = p_vals, T_vals = T_vals, 
              df = df)
  return(ret)
}
