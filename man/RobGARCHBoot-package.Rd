\name{RobGARCHBoot-package}
\alias{RobGARCHBoot-package}
\docType{package}
\title{
  Robust bootstrap forecast densities for GARCH models
}
\description{
  Bootstrap forecast densities for returns and volatilities using the robust residual-based bootstrap procedure of Trucíos et at. (2017). The package also includes the robust GARCH estimator of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2017).
}
\details{
This package provides a robust bootstrap procedure to obtain forecast densities for both return and volatilities in a GARCH context. The forecast densities are usefull to obtain forecast intervals as well as to estimate risk measures such as Value-at-Risk (VaR).
Additionally,  we also provide the robust GARCH estimator of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2017). 
This procedures showed good finite sample properties in both Monte Carlo experiments and empirial data. For a recent implementation of this procedure see Trucíos (2019).
}
\author{
Carlos Trucíos <ctrucios@gmail.com>
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
  
  
  Trucíos, Carlos. Forecasting Bitcoin risk measures: A robust approach. International Journal of Forecasting 35.3 (2019): 836-847.
}
\keyword{Bootstrap, Forecast Intervals, GARCH, Jump, Outliers, Robustness, Value-at-Risk (VaR), Volatility}
