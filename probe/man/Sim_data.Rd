\name{Sim_data}
\alias{Sim_data}
\docType{data}
\title{
Simulated high-dimensional data set for sparse linear regression
}
\description{
This dataset was simulated using a \eqn{100 \times 100} 2-dimensional setting described in the reference. The data contains 400 subjects with one outcome and 10,000 predictor variables. The test outcomes and predictor variables are contained in \code{Sim_data_test}.
}
\usage{data("Sim_data")}
\format{
  A data frame with 400 observations and the following objects:
  \describe{
    \item{\code{Y}}{Outcome variable of length \eqn{400}.}
    \item{\code{X}}{A \eqn{400 \times 10000} matrix of binary predictor variables.}
    \item{\code{signal}}{The locations of the non-zero regression coefficients.}
    \item{\code{beta_tr}}{The true values of all \eqn{10000} regression coefficients.}
    \item{\code{sigma2_tr}}{The true value of the residual variance.}
  }
}
\source{Simulated data.}
\examples{
data(Sim_data)
attach(Sim_data)
length(Y)
dim(X)
}
\keyword{datasets}
