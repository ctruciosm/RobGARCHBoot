\name{loglik_cDCC}
\alias{loglik_cDCC}
\title{
Loss function used in cDCC robust estimation.
}
\description{
Loss function used in cDCC robust estimation.
}
\usage{
loglik_cDCC(par, Qb, s, sigma)
}
\arguments{
  \item{par}{
teste
}
  \item{Qb}{
teste
}
  \item{s}{
teste
}
  \item{sigma}{
teste
}
}
\details{
This function is used in the robust estimation. We can use it to evaluate the value of the loss function using several values of the vector parameters (theta)
}
\value{
Returns the value of the loss function
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, C., Hotta, L. K., and Ruiz, E. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and value-at-risk. Journal of Statistical Computation and Simulation, 88.10 (2018): 1976-2000.
}
\author{
Carlos Trucíos
}
\examples{
\donttest{
# Estimating the parameters of the cDCC model in a robust way.
cDCC = Robust_cDCC(returns3)
param = cDCC[[1]]
Qbar = cDCC[[2]]
vol1 = fitted_Vol(param[1:3],returns3[,1])
vol2 = fitted_Vol(param[4:6],returns3[,2])
vol3 = fitted_Vol(param[7:9],returns3[,3])
e = matrix(c(returns3[,1]/vol1[1:nrow(returns3)],
             returns3[,2]/vol2[1:nrow(returns3)],
             returns3[,3]/vol3[1:nrow(returns3)]), ncol=3)
      
loglik_cDCC(param[10:11],Qbar,e, 0.8309765)      
}
}
\keyword{GARCH}
\keyword{Jump}
\keyword{Loss function}
\keyword{Outliers}
\keyword{Robustness}