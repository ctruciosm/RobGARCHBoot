#' @export
#' @import Rcpp
RobGARCHBoot <-
function(data, n.boot = 1000, n.ahead = 1){

# -----------------------
  coeff = ROBUSTGARCH(data)
  epps_c = scale(resBoot(coeff, data, coeff[1]/(1-coeff[2]-coeff[3]), 9),center = TRUE, scale = FALSE)

# -----------------------
	N = length(data); n.coeff = length(coeff)
	Coeff_b = matrix(NA, ncol = n.coeff, nrow = n.boot)                                                                                                    
	y_boot = matrix(NA, ncol = n.boot, nrow = N)                                                 
	yp = s2p = matrix(NA, ncol = n.ahead, nrow = n.boot)
		
		
# -----------------------
for (b in 1:n.boot) {
  e_b = rep(NA, n.ahead) 
  eboot = sample(epps_c, N, replace = TRUE)       
  
  # -----------------------
  y_boot[,b] = retBoot(coeff, coeff[2] + coeff[3], eboot,9)[[1]]
  coeffb =  ROBUSTGARCH(y_boot[,b])
  while(length(coeffb)!= n.coeff){
    eboot = sample(epps_c, N, replace = TRUE)                                      
    y_boot[,b] = retBoot(coeff,  coeff[2] + coeff[3],  eboot,9)[[1]]
    coeffb =  ROBUSTGARCH(y_boot[,b])
  }
  
  Coeff_b[b,] = coeffb
  
  # -----------------------
  s2_boot_B = sigma2Boot(Coeff_b[b,] , epps_c, Coeff_b[b,2] + Coeff_b[b,3], data, 9)

  
  # -----------------------
  e_b = sample(epps_c, n.ahead, replace = TRUE)
  forecast = foreBoot(Coeff_b[b,],  e_b, epps_c,  s2_boot_B , data,  n.ahead, 9)
   
  yp[b,1:n.ahead] = forecast[[1]][(N+1):(N+n.ahead)]
  s2p[b,1:n.ahead]= forecast[[2]][(N+1):(N+n.ahead)]

}
return(list(yp,s2p))
}
