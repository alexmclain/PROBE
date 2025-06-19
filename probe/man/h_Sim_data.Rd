\name{h_Sim_data}
\alias{h_Sim_data}
\docType{data}
\title{
Simulated high-dimensional data set for sparse linear regression with heterogeneous errorss
}
\description{
This dataset was simulated using a \eqn{20 \times 20} 2-dimensional setting described in the reference. The data contains 200 subjects with one outcome, 400 predictor variables, and 6 variables that are related to the variance of the error term. The data contains training and test data, along with the true values of the parameters used to generate the data..
}
\usage{data("h_Sim_data")}
\format{
  A data frame with 200 observations and the following objects:
  \describe{
    \item{\code{Y}}{Outcome variable of length \eqn{200}.}
    \item{\code{X}}{A \eqn{200 \times 400} matrix of binary predictor variables.}
    \item{\code{V}}{A \eqn{200 \times 7} design matrix to model variance of the residuals (includes an intercept).}
    \item{\code{beta_tr}}{The true values of all \eqn{400} regression coefficients.}
    \item{\code{omega_tr}}{The true values of all \eqn{7} variance regression coefficients.}
  }
}
\source{Simulated data.}
\examples{
data(h_Sim_data)
attach(h_sim_data)
length(Y)
dim(X)
}
\keyword{datasets}
