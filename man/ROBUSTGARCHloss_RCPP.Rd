\name{ROBUSTGARCHloss_RCPP}
\alias{ROBUSTGARCHloss_RCPP}
\title{
Loss function used in GARCH robust estimation.
}
\description{
Loss function used in GARCH (Generalized Autoregressive Conditional Heteroskedastic) robust estimation.
}
\usage{
ROBUSTGARCHloss_RCPP(theta, r, sigma2)
}
\arguments{
  \item{theta}{
Vector of robust estimated (or initial values) parameters obtained from ROBUSTGARCH function.
}
  \item{r}{
Vector of time series returns.
}
  \item{sigma2}{
robust squared volatility estimation (or inital value of squared volatility)
}
}
\details{
This functions is used in the robust estimation. We can use it to evaluate the value of the loss function using several values of the vector parameters (theta)
}
\value{
Returns the value of the loss function
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
# We can evaluate the loss function using the estimated parameters
ROBUSTGARCHloss_RCPP(param[2:3], returnsexample, param[1]/(1-param[2]-param[3]))
}

\keyword{GARCH}
\keyword{Jump}
\keyword{Loss function}
\keyword{Outliers}
\keyword{Robustness}