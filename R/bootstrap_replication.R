#' @noRd
#' @import Rcpp
bootstrap_replication <- function(data, epps_c, coeff, n.coeff, n.ahead, N){
  
  e_b = rep(NA, n.ahead) 
  eboot = sample(epps_c, N, replace = TRUE)       
  
  # -----------------------
  y_boot = retBoot(coeff, coeff[2] + coeff[3], eboot, 9)[[1]]
  coeffb =  ROBUSTGARCH(y_boot)
  while(length(coeffb)!= n.coeff){
    eboot = sample(epps_c, N, replace = TRUE)                                      
    y_boot = retBoot(coeff,  coeff[2] + coeff[3],  eboot,9)[[1]]
    coeffb =  ROBUSTGARCH(y_boot)
  }
  
  # -----------------------
  s2_boot_B = sigma2Boot(coeffb , epps_c, coeffb[2] + coeffb[3], data, 9)
  
  # -----------------------
  e_b = sample(epps_c, n.ahead, replace = TRUE)
  forecast = foreBoot(coeffb,  e_b, epps_c,  s2_boot_B , data,  n.ahead, 9)
  rs = c(forecast[[1]][(N+1):(N+n.ahead)],forecast[[2]][(N+1):(N+n.ahead)])
  
  return(rs)
}