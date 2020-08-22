\name{loglik_cDCC}
\alias{loglik_cDCC}
\title{
Likelihood function of the robust cDCC model
}
\description{
Loss function used in cDCC robust estimation.
}
\usage{
loglik_cDCC(par,Qb,s, sigma)
}
\arguments{
  \item{par}{
Vector of cDCC parameters
}
  \item{Qb}{
Qbar matrix in the cDCC equation
}
\item{s}{
matrix of residuals coming from the univariate robust GARCH estimation.
}
  \item{sigma}{
parameter used in Equation (2.8) of Boudt et al. (2013) 
}
}
\details{
This function is used in the robust cDCC estimation. 
}
\value{
Returns the value of the loss function
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and value-at-risk. Journal of Statistical Computation and Simulation 88.10 (2018): 1976-2000.
}
\author{
Carlos Trucíos
}
\examples{
\donttest{
loglik_cDCC(par,Qb, s, sigma)
}
}

\keyword{cDCC}
\keyword{Jump}
\keyword{Loss function}
\keyword{Outliers}
\keyword{Robustness}