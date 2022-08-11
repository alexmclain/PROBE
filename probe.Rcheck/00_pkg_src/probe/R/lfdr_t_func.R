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
  y[y<=0] <- 1e-10
  lfdr <- pi0 * dt(T,df=df_val)/y
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
  res <- list(lfdr=lfdr_out,f=f)
  return(res)
}
