#' @details Examples for applying PROBE to sparse high-dimensional linear regression are given
#' for one-at-a-time \code{\link{probe_one}} or all-at-once \code{\link{probe}} type optimization.
#'
#' @references \itemize{ \item McLain, AC, A Zgodic, H Bondell (2025). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. \textit{Computational Statistics and Data Analysis} 207, 108146.
#' \item Zgodic, A., Bai, R., Zhang, J., Wang, Y., Rorden, C., & McLain, A. (2023). Quantifying predictive uncertainty of aphasia severity in stroke patients with sparse heteroscedastic Bayesian high-dimensional regression. arXiv preprint arXiv:2309.08783.}
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