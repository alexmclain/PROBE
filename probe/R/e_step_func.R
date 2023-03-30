#' @title Function for fitting the empirical Bayes portion of the E-step
#' @description A wrapper function estimating posterior expectations of the \eqn{\gamma} variables using an empirical Bayesian technqiue. 
#' @param beta_t Expectation of the posterior mean (assuming \eqn{\gamma=1})
#' @param beta_var Current posterior variance (assuming \eqn{\gamma=1})
#' @param df Degrees of freedom for the t-distribution (use to calculate p-values).
#' @param adj Bandwidth multiplier to Silverman's `rule of thumb' for calculating the marginal density of the test-statistics (default = 5).
#' @param lambda Value of the \eqn{\lambda} parameter for estimating the proportion of null hypothesis using Storey et al. (2004) (default = 0.1).
#' @param monotone Logical - Should the estimated marginal density of the test-statistics be monotone non-increasing from zero (default = TRUE).
#' @return A list including 
#' 
#' \code{delta} estimated posterior expectations of the \eqn{\gamma}.
#' 
#' \code{pi0} estimated proportion of null hypothesis
#' @references 
#' Storey, J. D., Taylor, J. E., and Siegmund, D. (2004), “Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: A unified approach,” J. R. Stat. Soc. Ser. B. Stat. Methodol., 66, 187–205.
#' McLain, A. C., Zgodic, A., & Bondell, H. (2022). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2209.08139.
#' @examples
#' #not run
#' #mod <- e_step_func(beta_t, beta_var, df, adj = 5, lambda = 0.1, monotone = TRUE)

e_step_func <- function(beta_t, beta_var, df, adj = 5, lambda = 0.1, monotone=TRUE){
  
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
  EB_res <- lfdr_t_func(T = T_vals, pi0 = pi_hat, trunc=TRUE, 
                         monotone = monotone, adj = adj, 
                         df_val = df)
  gamma <- 1-EB_res$lfdr
  adj_warning <- EB_res$adj_warning
  
  if(sum(gamma)==0 & adj_warning == 1){
    gamma <- 1-lfdr_t_func(T = T_vals, pi0 = pi_hat, trunc=TRUE, 
                           monotone = TRUE, adj = 2*adj, 
                           df_val = df)$lfdr
    if(sum(gamma)==0){
      gamma <- 1-lfdr_t_func(T = T_vals, pi0 = pi_hat, trunc=TRUE, 
                             monotone = FALSE, adj = 2*adj, 
                             df_val = df)$lfdr
    }
  }
  
  ret <- list(beta_tilde = beta_t, beta_tilde_var = beta_var, gamma = gamma, 
              lfdr = 1 - gamma, pi0 = pi_hat, p_vals = p_vals, T_vals = T_vals, 
              df = df, adj_warning = adj_warning)
  return(ret)
}

pi0_func <- function (p, lambda = 0.1) 
{
  rm_na <- !is.na(p)
  p <- p[rm_na]
  m <- length(p)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  }
  
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  
  return(list(pi0 = pi0))
}


lfdr_t_func <- function (T, pi0 = NULL, trunc = TRUE, monotone = TRUE, adj=3, 
                         df_val = NULL, one.sided = FALSE, ...) 
{
  lfdr_out <- T
  rm_na <- !is.na(T)
  T <- T[rm_na]
  if (is.null(pi0)) {
    stop("pi0 must be given.")
  }
  n <- length(T)
  myd <- density(T, adjust = adj)
  mys <- smooth.spline(x = myd$x, y = myd$y)
  y <- predict(mys, T)$y
  f0 <- dt(T,df=df_val)
  max_d <- c(max(y), predict(mys, 0)$y)
  max_n <- c(max(f0), dt(0,df=df_val))
  y[y<=0] <- 1e-10
  lfdr <- pi0 * f0/y
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    T_neg <- T[T<0]
    lfdr_neg <- lfdr[T<0]
    T_pos <- T[T>=0]
    lfdr_pos <- lfdr[T>=0]
    
    o <- order(T_neg, decreasing = FALSE)
    ro <- order(o)
    lfdr[T<0] <- cummax(lfdr_neg[o])[ro]
    
    o <- order(T_pos, decreasing = TRUE)
    ro <- order(o)
    lfdr[T>=0] <- cummax(lfdr_pos[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  f <- list(x=sort(T),y=lfdr[order(T)])
  
  adj_warning = 0
  if(any(max_d > max_n)){
    adj_warning = 1
  }
  res <- list(lfdr = lfdr_out, f = f, 
              adj_warning = adj_warning)
  return(res)
}
