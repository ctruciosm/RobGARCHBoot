#include<RcppArmadillo.h>
#include<Rmath.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
//' @export
//' @noRd
// [[Rcpp::export]]
SEXP ROBUSTGARCHloss_RCPP(NumericVector theta, NumericVector r, double sigma2){
  int n = r.size(), k = 3;
  NumericVector h(n), J(n), auxr(n), y(n);
  Rcpp::Function mean("mean");
  h[0] = sigma2;
  J[0] = r[0]/sqrt(h[0]);
  
  for(int i=1; i<n; i++){
    if(std::abs(J[i-1])<k){
      h[i]= sigma2*(1-theta[0]-theta[1])+ theta[0]*pow(r[i-1],2)+ theta[1]*h[i-1];
    } else{
      h[i]= sigma2*(1-theta[0]-theta[1])+ theta[0]*1.005018*h[i-1]+ theta[1]*h[i-1];
    }
    J[i] = r[i]/sqrt(h[i]);
  }
  auxr = ifelse(r==0,r+0.00001,r);
  y = log(pow(auxr,2)/h);
  return(wrap(mean(-y + 0.8260*5*log(1+exp(y)/2))));
}

//' @noRd
//' @useDynLib RobGARCHBoot
// [[Rcpp::export]]
SEXP grid_RCPP(NumericVector y, double sigmaR){
  NumericVector coeff(2),vi(2);
  double alfa1, beta1;
  double alfa1min = 0.005,  alfa1max = 0.2, beta1min = 0.65,  beta1max = 0.98,  nalfa1 = 5,  nbeta1 = 5;
  double ml = 100000000, nml;
  double lmalfa1 = (alfa1max-alfa1min)/nalfa1;
  double lmbeta1 = (beta1max-beta1min)/nbeta1;
  Rcpp::Function ROBUSTGARCHloss_RCPP("ROBUSTGARCHloss_RCPP");
  
  for(int nj=0; nj<nalfa1; nj++){
    for(int nk=0; nk<nbeta1; nk++){
      alfa1 = alfa1min+nj*lmalfa1;
      beta1 = beta1min+nk*lmbeta1; 
      
      if(alfa1+beta1<0.999){
        coeff[0] = alfa1;
        coeff[1] = beta1;
        nml=Rcpp::as<double>(ROBUSTGARCHloss_RCPP(coeff,y,sigmaR));
        if (nml<ml){
          vi[0] = coeff[0];
          vi[1] = coeff[1];
          ml=nml;
        }
      }
    }
  }
  return(vi);
}

//' @noRd
//' @useDynLib RobGARCHBoot
// [[Rcpp::export]]
SEXP resBoot(NumericVector coeff, NumericVector r, double S, double k){
  int n = r.size();
  NumericVector h(n), e(n), aux(n-1);
  h[0] = S;
  e[0] = r[0]/sqrt(h[0]);
  for(int i=1; i<n; i++){
    aux[i-1] = pow(r[i-1],2)/h[i-1];
    if(aux[i-1]<k){
      h[i] = coeff[0] + coeff[1]*pow(r[i-1],2) + coeff[2]*h[i-1];
      e[i] = r[i]/sqrt(h[i]);
    } else {
      h[i] = coeff[0] + (coeff[1]*1.005018+coeff[2])*h[i-1];
      e[i] = r[i]/sqrt(h[i]);
    }
  }
  return(e);
}

//' @noRd
//' @useDynLib RobGARCHBoot
// [[Rcpp::export]]
SEXP retBoot(NumericVector coeff, double S, NumericVector e, double k){
  int n = e.size();
  double suncond;
  Rcpp::Function sample("sample");
  NumericVector e_s(1);
  NumericVector r(n), h(n), aux(n);
  suncond = coeff[0]/(1-S);
  
  r[0] = e[0]*sqrt(suncond);
  h[0] = suncond;
  
  for(int i=1; i<n; i++){
    aux[i-1] = pow(r[i-1],2)/h[i-1];
    if(aux[i-1]<k){
      h[i] = coeff[0] + coeff[1]*pow(r[i-1],2) + coeff[2]*h[i-1];
      r[i] = e[i]*sqrt(h[i]);
    } else {
      e_s[0]= Rcpp::as<double>(sample(e,1));
      h[i] = coeff[0] + (coeff[1]*pow(e_s[0],2)+coeff[2])*h[i-1];
      r[i] = e[i]*sqrt(h[i]);
    }
  }
  return Rcpp::List::create(r,h);
}

//' @noRd
//' @useDynLib RobGARCHBoot
// [[Rcpp::export]]
SEXP sigma2Boot(NumericVector coeff, NumericVector e, double S, NumericVector r, double k){
  Rcpp::Function sample("sample");
  NumericVector e_s(1);
  int n = r.size();
  double suncond;
  NumericVector h(n), aux(n-1);
  suncond = coeff[0]/(1-S);
  h[0] = suncond;
  for(int i=1; i<n; i++){
    aux[i-1] = pow(r[i-1],2)/h[i-1];
    if(aux[i-1]<=k){
      h[i] = coeff[0] + coeff[1]*pow(r[i-1],2) + coeff[2]*h[i-1];
    } else {
      e_s[0]= Rcpp::as<double>(sample(e,1));
      h[i] = coeff[0] + (coeff[1]*pow(e_s[0],2)+coeff[2])*h[i-1];
    }
  }
  return(h);
}

//' @noRd
//' @useDynLib RobGARCHBoot
// [[Rcpp::export]]
SEXP foreBoot(NumericVector coeff, NumericVector e, NumericVector e2, NumericVector h, NumericVector r, int ahead, double k){
  int n = r.size();
  NumericVector hp(n+ahead), rp(n+ahead), aux(ahead);
  Rcpp::Function sample("sample");
  NumericVector e_s(1);
  
  for(int j = 0; j<n; j++){
    hp[j] = h[j];
    rp[j] = r[j];
  }
  for(int i = 0; i<ahead; i++){
    aux[i] = pow(rp[n-1+i],2)/hp[n-1+i];
    if(aux[i]<=k){
      hp[n+i] = coeff[0] + coeff[1]*pow(rp[n-1+i],2) + coeff[2]*hp[n-1+i];
      rp[n+i] = e[i]*sqrt(hp[n+i]);
    } else {
      e_s[0]= Rcpp::as<double>(sample(e2,1));
      hp[n+i] = coeff[0] + (coeff[1]*pow(e_s[0],2)+coeff[2])*hp[n-1+i];
      rp[n+i] = e[i]*sqrt(hp[n+i]);
    }
  }
  return Rcpp::List::create(rp,hp,aux);
}



// [[Rcpp::depends(RcppArmadillo)]]
//' @noRd
// [[Rcpp::export]]
SEXP gridcDCC(arma::mat Qb,arma::mat s, double sigma){
  NumericVector coeff(2),vi(2);
  double alfa1, beta1;
  double alfa1min = 0.01,  alfa1max = 0.3, beta1min = 0.65,  beta1max = 0.98,  nalfa1 = 5,  nbeta1 = 5;
  double ml = 100000000, nml;
  double lmalfa1 = (alfa1max-alfa1min)/nalfa1;
  double lmbeta1 = (beta1max-beta1min)/nbeta1;
  Rcpp::Function loglik_cDCC("loglik_cDCC");
  
  for(int nj=0; nj<nalfa1; nj++){
    for(int nk=0; nk<nbeta1; nk++){
      alfa1 = alfa1min+nj*lmalfa1;
      beta1 = beta1min+nk*lmbeta1; 
      
      if(alfa1+beta1<0.999){
        coeff[0] = alfa1;
        coeff[1] = beta1;
        nml=Rcpp::as<double>(loglik_cDCC(coeff,Qb,s, sigma));
        if (nml<ml){
          vi[0] = coeff[0];
          vi[1] = coeff[1];
          ml=nml;
        }
      }
    }
  }
  return(vi);
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @export
//' @noRd
// [[Rcpp::export]]
SEXP loglik_cDCC(arma::vec par,arma::mat Qb,arma::mat s, double sigma){
  int T = s.n_rows;
  int K = s.n_cols;
  Rcpp::Function pchisq("pchisq");
  Rcpp::Function qchisq("qchisq");
  double CN = K/(K*as<double>(pchisq(qchisq(0.9973,K),K+2)) + as<double>(qchisq(0.9973,K))*(1-0.9973));
  arma::vec lR1(T), lR2(T), lR3(T), d(T);
  arma::mat  R(K,K), Qt(K,K),Pt(K,K), iPt(K,K);
  R.zeros();  Qt.zeros();  Pt.zeros();
  R = Qb;
  Qt = Qb;
  lR1[0] = log(det(R));
  d[0] = arma::conv_to<double>::from(s.row(0)*R.i()*s.row(0).t());
  lR2[0] = (K+4)*sigma*log(1+d[0]/2);
  lR3[0] =lR2[0]+lR1[0];
  Pt = sqrt(diagmat(Qt));
  
  for(int t = 1; t<T; t++){
    if (1<(as<double>(qchisq(0.9973,K))/d[t-1])){ 
      Qt = (1-par[0]-par[1])*Qb+ par[0]*CN*Pt*s.row(t-1).t()*s.row(t-1)*Pt +par[1]*Qt;
    } else {
      Qt = (1-par[0]-par[1])*Qb + par[0]*CN*K/d[t-1]*Pt*s.row(t-1).t()*s.row(t-1)*Pt +par[1]*Qt;
    }
    Pt = sqrt(diagmat(Qt));
    iPt = Pt.i();
    R =  iPt*Qt*iPt;
    lR1[t] = log(det(R));
    d[t] = arma::conv_to<double>::from(s.row(t)*R.i()*s.row(t).t());
    lR2[t] = (K+4)*sigma*log(1+d[t]/2); 
    lR3[t] = lR2[t]+lR1[t];
  }
  return Rcpp::wrap(mean(lR3));
}
