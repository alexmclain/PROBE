\name{Sim_data_cov}
\alias{Sim_data_cov}
\docType{data}
\title{
Simulated high-dimensional data set for sparse linear regression with non-sparse covariates.
}
\description{
This dataset was simulated using a \eqn{100 \times 100} 2-dimensional setting described in the reference only two covariates are added. The data contains 400 subjects with one outcome, 10,000 predictor variables which are to be subjected to the sparsity assumption, and 2 covariates which are not to be subjected to the sparsity assumption. 
}
\usage{data("Sim_data_cov")}
\format{
  A data frame with 400 observations and the following objects:
  \describe{
    \item{\code{Y}}{Outcome variable of length \eqn{400}.}
    \item{\code{X}}{A dataframe of a continuous (\code{Cont_cov}) and binary (\code{Binary_cov}) covariate.}
    \item{\code{Z}}{A \eqn{400 \times 10000} matrix of binary predictor variables.}
    \item{\code{beta_tr}}{The true values of all \eqn{10000} regression coefficients.}
    \item{\code{beta_X_tr}}{The true values of the intercept, \code{Cont_cov}, and \code{Binary_cov}.}
    \item{\code{signal}}{The locations of the non-zero regression coefficients.}
  }
}
\references{
To come
}
\examples{
data(Sim_data_cov)
attach(Sim_data_cov)
length(Y)
summary(X)
dim(Z)
}
\keyword{datasets}
