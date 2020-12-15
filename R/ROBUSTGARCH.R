#' @noRd
#' @importFrom stats median
medianB = function(r){
  k = 30
  n = length(r)
  med = c()
  mad = c()
  for(i in 1:n){
    med[i] = ifelse(i<k/2,median(r[1:(k+1)]),ifelse(i>n-k/2,median(r[(n-k):n]),median(r[(i-k/2):(i+k/2)])))
    mad[i] = ifelse(i<k/2,median(abs(r[1:(k+1)] - med[i])),ifelse(i>n-k/2,median(abs(r[(n-k):n] - med[i])),median(abs(r[(i-k/2):(i+k/2)] - med[i]))))
  }
  return(list(med,mad))
}

#' @export
fitted_Vol = function(theta,r){
  n= length(r)+1
  h= c()
  k = 3
  h[1]= theta[1]/(1-theta[2]-theta[3])
  for (t in 2:n){
    if(abs(r[t-1]/sqrt(h[t-1]))<k){
      h[t]= theta[1]+ theta[2]*r[t-1]^2+ theta[3]*h[t-1]
    } else{
      h[t]= theta[1]+ theta[2]*1.005018*h[t-1]+ theta[3]*h[t-1]
    }
  }
  return(sqrt(h))
}

#' @export
#' @import Rcpp
#' @importFrom stats constrOptim
ROBUSTGARCH = function(y){
  AUX= medianB(y)
  Med = AUX[[1]]
  MAD = AUX[[2]]
  I = (y-Med)^2/(1.486*MAD)^2<= 3.841459
  mu_R = sum(y*I)/sum(I)
  J = (y-mu_R)^2/(1.486*MAD)^2<=3.841459
  sigma2R = 1.318 * sum((y-mu_R)^2*J)/sum(J)
  ini = grid_RCPP(y-mu_R, sigma2R)
  ra <- matrix(c(1,0,0,1,-1,-1),ncol=2,byrow=TRUE)
  rb <- c(0.00001, 0.00001,-0.9999)
  param <- constrOptim(theta = ini, f = ROBUSTGARCHloss_RCPP, grad = NULL, ui = ra, ci = rb, sigma2 = sigma2R, r = y, outer.iterations = 400, outer.eps = 1e-07)$par
  param = c(sigma2R*(1-param[1]-param[2]),param[1],param[2])
  return(param)
}


#' @noRd
#' @importFrom stats cor qchisq pchisq
Q_bar = function(s, k = 30, gamma = 0.95){
  n = nrow(s)
  p = ncol(s)
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
  n = nrow(r)
  p = ncol(r)
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
  
  ra = matrix(c(1,0,0,1,-1,-1),ncol=2,byrow=TRUE)
  rb = c(0.00001, 0.00001,-0.9999)
  coef_cDCC = constrOptim(theta = parini, f = loglik_cDCC, grad = NULL, ui = ra, ci = rb, Qb = Qbarra,s = e, sigma = sigma_)$par
  
  coeff = c(as.vector(coef),coef_cDCC)
  return(list(coeff,Qbarra))
}


#' @noRd
#' @import Rcpp
fitted_cDCC = function(r, Qbar, params){
  Dim = dim(r)
  n = Dim[1]
  p = Dim[2]
  H = list()
  vol = matrix(0,ncol = p, nrow = n+1)
  e = matrix(0, ncol = p, nrow = n)
  for (i in 1:p){
    coef = params[(3*(i-1)+1):(3*i)]
    vol[,i] = fitted_Vol(coef, r[,i])
    e[,i] = r[,i]/vol[1:n,i]
  }
  
  dccpar = params[(3*p+1):(3*p+2)]
  R = cor_cDCC(dccpar,Qbar,e)
  
  for(j in 1:(n+1)){
    H[[j]] = diag(vol[j,])%*%R[,,j]%*%diag(vol[j,])
  }
  
  return(H)
}
