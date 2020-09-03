\name{fitted_cDCC}
\alias{fitted_cDCC}
\title{
Estimated conditional covariance matrix.
}
\description{
Using the robust procedure of Boudt et al. (2013) with the modifications introduced by Trucíos et at. (2018), we obtain the estimated conditional covariance matrix.
}
\usage{
fitted_cDCC(r, Qbar, params)
}
\arguments{
  \item{r}{
matrix of time series returns
}
\item{Qbar}{
Qbar matrix in the cDCC equation
}
\item{params}{
Vector of parameters
}
}
\details{
More details can be found in Boudt et al. (2013) and Trucíos et at. (2018).
}
\value{
The function returns the estimated conditional covariance matrix from 1 to T+1.
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
cDCC = Robust_cDCC(r = cryptoreturn)
H = fitted_cDCC(r = cryptoreturn, Qbar = cDCC[[2]], params = cDCC[[1]])
}
}

\keyword{cDCC}
\keyword{Covariance}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}