#' @noRd
#' @importFrom stats cor qchisq pchisq
Q_bar = function(s, k = 30, gamma = 0.95){
  Dimension = dim(s)
  n = Dimension[1]
  p = Dimension[2]
  AUXL = 20
  Lt = rep(0,n)
  Su = matrix(0,p,p)
  for(i in 1:n){
    if (i<=k/2){
      Ct = cor(s[1:(k+1),],method="spearman")
    } else {
      if (i>n-k/2){
        Ct = cor(s[(n-k):n,],method="spearman")
      } else {
        Ct = cor(s[(i-k/2):(i+k/2),],method="spearman")
      }
    }
    SCt = 2*sin(1/6*pi*Ct)
    AUXL = try(t(s[i,])%*%solve(SCt)%*%s[i,],silent = TRUE)
    Lt[i] = ifelse( AUXL <= qchisq(gamma,p),1,0)       
    Su = Su + s[i,]%*%t(s[i,])*Lt[i]
  }
  CN = (p/(p*pchisq(qchisq(gamma,p),p+2) + (1-gamma)*qchisq(gamma,p)))
  RC = CN*Su/sum(Lt)  
  D = solve(sqrt(diag(diag(RC))))
  R = D%*%RC%*%D
  return(R)
}


#' @export
#' @import Rcpp
#' @importFrom stats constrOptim pchisq qchisq integrate dchisq
Robust_cDCC = function(r){
  Dim = dim(r)
  n = Dim[1]
  p = Dim[2]
  coef = matrix(0, ncol = p, nrow = 3)
  vol = matrix(0,ncol = p, nrow = n+1)
  e = matrix(0, ncol = p, nrow = n)
  for (i in 1:p){
    coef[,i] = ROBUSTGARCH(r[,i])
    vol[,i] = fitted_Vol(coef[,i], r[,i])
    e[,i] = r[,i]/vol[1:n,i]
  }
  
  integrand = function(x) { (p+4)*x/(2+x)*dchisq(x,p) }
  sigma_ = p/integrate(integrand,0,Inf)$value
  
  Qbarra = Q_bar(e)
  parini = gridcDCC(Qbarra,e, sigma_)
  
  loglik_cDCC(parini, Qbarra, e, sigma_)
  
  ra = matrix(c(1,0,0,1,-1,-1),ncol=2,byrow=TRUE)
  rb = c(0.00001, 0.00001,-0.9999)
  coef_cDCC = constrOptim(theta = parini, f = loglik_cDCC, grad = NULL, ui = ra, ci = rb, Qb = Qbarra,s = e, sigma = sigma_)$par
  
  coeff = c(as.vector(coef),coef_cDCC)
  return(list(coeff,Qbarra))
}