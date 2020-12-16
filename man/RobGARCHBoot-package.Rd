\name{RobGARCHBoot-package}
\alias{RobGARCHBoot-package}
\docType{package}
\title{
  Robust Bootstrap Forecast Densities for GARCH  Models
}
\description{
Bootstrap forecast densities for returns and volatilities using the robust residual-based bootstrap procedure of Trucíos et at. (2017). The package also includes the robust GARCH (Generalized Autoregressive Conditional Heteroskedastic) estimator of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2017). The robust cDCC estimator used in Trucíos et at. (2018) is also implemented
}
\details{
This package provides a robust bootstrap procedure to obtain forecast densities for both return and volatilities in a GARCH context. The forecast densities are useful to obtain forecast intervals as well as to estimate risk measures such as Value-at-Risk (VaR) and Expected Shortfall (ES).
We also provide the robust GARCH estimator of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2017). 
This procedure has shown good finite sample properties in both Monte Carlo experiments and empirical data. See;  Trucíos et al. (2017), Trucíos (2019) and Trucíos et al. (2020) for recent implementations.

}
\author{
Carlos Trucíos <ctrucios@gmail.com>
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
  
  Trucíos, C., Hotta, L. K., and Ruiz, E. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and value-at-risk. Journal of Statistical Computation and Simulation, 88.10 (2018): 1976-2000.
  
  Trucíos, Carlos. Forecasting Bitcoin risk measures: A robust approach. International Journal of Forecasting 35.3 (2019): 836-847.
  
  Trucíos, Carlos, Aviral K. Tiwari, and Faisal Alqahtani. Value-at-risk and expected shortfall in cryptocurrencies' portfolio: a vine copula–based approach. Applied Economics 52.24 (2020): 2580-2593.
}
\keyword{Bootstrap}
\keyword{DCC}
\keyword{Forecast Intervals}
\keyword{GARCH}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}
\keyword{Value-at-Risk (VaR)}
\keyword{Volatility}
