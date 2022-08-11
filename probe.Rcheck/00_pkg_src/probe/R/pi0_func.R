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