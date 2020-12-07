#' @export
#' @import Rcpp foreach doParallel doRNG
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
RobGARCHBootParallel <-
function(data, n.boot = 1000, n.ahead = 1){

# -----------------------
  coeff = ROBUSTGARCH(data)
  epps_c = scale(resBoot(coeff, data, coeff[1]/(1-coeff[2]-coeff[3]), 9),center = TRUE, scale = FALSE)

# -----------------------
	N = length(data)
  n.coeff = length(coeff)

# -----------------------
	
 cl <- makeCluster(detectCores(logical = FALSE)-1, setup_strategy = "sequential")
 registerDoParallel(cl)	
 rsboot = foreach(b = 1:n.boot, .combine = rbind)	%dorng% {
   bootstrap_replication(data, epps_c, coeff, n.coeff, n.ahead, N)
 }
 stopCluster(cl)
 yp  =  matrix(rsboot[,1:n.ahead], ncol = n.ahead)
 s2p =  matrix(rsboot[,(n.ahead+1):(n.ahead*2)], ncol = n.ahead)

 return(list(yp,s2p))
}
