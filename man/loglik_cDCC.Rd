\name{loglik_cDCC}
\alias{loglik_cDCC}
\title{
Loss function used in cDCC robust estimation.
}
\description{
Loss function used in cDCC robust estimation.
}
\usage{
loglik_cDCC(par, Qb, s, sigma)
}
\arguments{
  \item{par}{
teste
}
  \item{Qb}{
teste
}
  \item{s}{
teste
}
  \item{sigma}{
teste
}
}
\details{
This function is used in the robust estimation. We can use it to evaluate the value of the loss function using several values of the vector parameters (theta)
}
\value{
Returns the value of the loss function
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, C., Hotta, L. K., & Ruiz, E. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and value-at-risk. Journal of Statistical Computation and Simulation, 88.10 (2018): 1976-2000.
}
\author{
Carlos Trucíos
}
\keyword{GARCH}
\keyword{Jump}
\keyword{Loss function}
\keyword{Outliers}
\keyword{Robustness}