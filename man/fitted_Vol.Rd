\name{fitted_Vol}
\alias{fitted_Vol}
\title{
Estimated Volatility
}
\description{
Using the robust estimated parameters of Boudt et al. (2013) and Trucíos et at. (2017), we obtain the estimated volatility.
}
\usage{
fitted_Vol(theta,r)
}
\arguments{
  \item{theta}{
Vector of robust estimated parameters obtained from ROBUSTGARCH function.
}
  \item{r}{
Vector of time series returns.
}
}
\details{
More details can be found in Boudt et al. (2013) and Trucíos et at. (2017).
}
\value{
The function returns the estimated volatility from 1 to T+1.
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
}
\author{
Carlos Trucíos
}
\examples{
# Using the Bitcoin daily returns, we estimate the parameter of the GARCH model in a robust way
param = ROBUSTGARCH(returnsexample)
# With the estimated parameters, we estimate the volatiltiy in a robust way
vol = fitted_Vol(param, returnsexample)
}

\keyword{GARCH}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}