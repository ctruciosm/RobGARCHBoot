\name{Robust_cDCC}
\alias{Robust_cDCC}
\title{
Robust cDCC Estimator
}
\description{
Robust cDCC estimator as in Trucíos et at. (2018), which is a slight modification of the procedure of Boud et al. (2013).
}
\usage{
Robust_cDCC(r)
}
\arguments{
  \item{r}{
Matrix of time series returns.
}
}
\details{
More details can be found in Boudt et al. (2013) and Trucíos et at. (2018).
}
\value{
The function returns the estimated parameters and the Qbar matrix
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
# Estimating the parameters of the cDCC model in a robust way.
cDCC = Robust_cDCC(returns3)
#parameters
cDCC[[1]]
#Qbar
cDCC[[2]]
}
}

\keyword{DCC}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}