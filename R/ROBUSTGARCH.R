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



