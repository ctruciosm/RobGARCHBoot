\name{ROBUSTGARCH}
\alias{ROBUSTGARCH}
\title{
Robust GARCH estimator
}
\description{
Robust GARCH estimator of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2017).
}
\usage{
ROBUSTGARCH(y)
}
\arguments{
  \item{y}{
Vector of time series returns.
}
}
\details{
More details can be found in Boudt et al. (2013) and Trucíos et at. (2017).
}
\value{
The funcion returns the estimated parameters.
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
}
\author{
Carlos Trucíos
}
\examples{
param = ROBUSTGARCH(returnsexample*100)
param
}

\keyword{GARCH, Jump, Outliers, Robustness}