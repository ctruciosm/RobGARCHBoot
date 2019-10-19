\name{RobGARCHBoot}
\alias{RobGARCHBoot}
\title{
Robust GARCH Estimator
}
\description{
Robust Bootstrap procedure of Trucíos et al. (2017)
}
\usage{
RobGARCHBoot(data, n.boot = 1000, n.ahead = 1)
}
\arguments{
  \item{data}{
Vector of time series returns.
}
  \item{n.boot}{
Number of bootsrap replications. By default n.boot = 1000
}
  \item{n.ahead}{
Numbers of steps-ahead. By default n.ahead = 1
}

}
\details{
More details can be found in Trucíos et at. (2017), Hotta and Trucíos (2018) and Trucíos (2019).
}
\value{
The funcion returns two lists with the empirical H-steps-ahead bootstrap densities for returns and squared volatilities.
}
\references{
  Hotta, Luiz Koodi, and Carlos Trucíos. Inference in (M)GARCH models in the presence of additive outliers: Specification, estimation, and prediction. Advances in Mathematics and Applications. Springer, Cham, 2018. 179-202.

  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
  
  Trucíos, Carlos. Forecasting Bitcoin risk measures: A robust approach. International Journal of Forecasting 35.3 (2019): 836-847.
}
\author{
Carlos Trucíos
}
\examples{
# boot = RobGARCHBoot(returnsexample, n.boot = 1000, n.ahead = 1)
# One-step-ahead VaR 1\%
# VaR1 = quantile(boot[[1]], prob = 0.01)
# Forecast Intervals 95\%: Returns
# quantile(boot[[1]], prob = c(0.025, 0.975))
# Forecast Intervals 95\%: Squared volatilities
# quantile(boot[[2]], prob = c(0.025, 0.975))
}
\keyword{Bootstrap}
\keyword{GARCH}
\keyword{Forecast Intervals}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}