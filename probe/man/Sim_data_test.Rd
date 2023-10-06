\name{Sim_data_test}
\alias{Sim_data_test}
\docType{data}
\title{
Simulated high-dimensional test data set for sparse linear regression
}
\description{
A test set of outcomes and predictor variables to be used with \code{Sim_data}.
}
\usage{data("Sim_data_test")}
\format{
  A data frame with 400 observations and the following objects:
  \describe{
    \item{\code{Y_test}}{Outcome variable of length \eqn{400} for test set.}
    \item{\code{Z_test}}{A \eqn{400 \times 10000} matrix of binary predictor variables  for test set.}
  }
}
\source{Simulated data.}
\examples{
data(Sim_data_test)
attach(Sim_data_test)
length(Y_test)
dim(X_test)
}
\keyword{datasets}
