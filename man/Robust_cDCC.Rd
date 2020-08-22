\name{Robust_cDCC}
\alias{Robust_cDCC}
\title{
Robust cDCC Estimator
}
\description{
Robust cDCC (Dynamic Conditional Correlation) estimator of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2018).
}
\usage{
Robust_cDCC(r)
}
\arguments{
  \item{r}{
Nxp matrix of time series returns where N stands for the sample size and p stands for the number of time series.
}
}
\details{
More details can be found in Boudt et al. (2013) and Trucíos et at. (2018).
}
\value{
The function returns the estimated parameters.
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
Robust_cDCC(y)
}
}

\keyword{cDCC}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}