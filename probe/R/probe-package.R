#' @details Examples for applying PROBE to sparse high-dimensional linear regression are given
#' for one-at-a-time \code{\link{probe_one}} or all-at-once \code{\link{probe}} type optimization.
#'
#' @references \itemize{ \item McLain, A. C., Zgodic, A., & Bondell, H. (2022). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2209.08139.}
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom glmnet cv.glmnet
#' @importFrom graphics axis legend points mtext par
#' @importFrom stats coef density  dt influence lm p.adjust pchisq predict pt qnorm quantile smooth.spline var vcov
#' @aliases probe-package
#' @useDynLib probe, .registration = TRUE
"_PACKAGE"


NULL