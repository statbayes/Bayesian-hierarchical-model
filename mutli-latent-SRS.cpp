//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include<RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;


//=========================================================================
// get upper/lower triangular without diagonal (get correlations vector)
//=========================================================================

// [[Rcpp::export]]

colvec corvec(mat V){
  int l = V.n_rows;
  colvec x(1);
  
  for(int j=0; j<l-1; j++){
    
    colvec temp = trans(V(j, span(j+1,l-1)));
    x = join_cols(x,temp);
    
  }
  
  x.shed_row(0);
  
  return x;
}

// [[Rcpp::export]]
colvec varvec(mat V){
  int l = V.n_rows;
  colvec x(1);
  
  for(int j=0; j<l; j++){
    
    colvec temp = trans(V(j, span(j,l-1)));
    x = join_cols(x,temp);
   
  }
  
  x.shed_row(0);
  
  return x;
}



// [[Rcpp::export]]
mat mvrnorm(int n, colvec mu, mat sigma) {
   int ncols = sigma.n_cols;
   mat Y = randn(n, ncols);
   return (repmat(mu, 1, n).t() + Y * chol(sigma));
}



//[[Rcpp::export]]
mat rwish(double nu, mat S) {
  //Rcpp::NumericVector nu(NU);
  //arma::mat S = as<arma::mat>(s);
  mat CC = chol(S);
  int n = S.n_cols;
  mat x = randn(n, n);
  
  for (int i = 0; i < n; i++) {
    x.diag()[i] = sqrt(R::rchisq(nu-i));
  }
  x = trimatu(x);
  x = x * CC;
  x = trans(x) * x;
  return x;
}

//[[Rcpp::export]]
mat riwish(double nu, mat S){
  int n = S.n_cols;
  mat X(n,n);
  X = inv(rwish(nu,inv(S)));
  return X;
}



//=============================================================
// check positive definite matrix
//=============================================================


// [[Rcpp::export]]
bool ispd(mat X, double tol=1e-8){
  colvec eigenval = eig_sym(X);
  
  int n = X.n_rows;
  //double maxeig = max(eigenval);
  //if(!all(eigenval>=-tol*abs(maxeig))){
  //  return false;
  //}
  colvec abseig = abs(eigenval);
  
  for(int i=0;i<n;i++){
    if(abseig(i) < tol){
      eigenval(i) = 0;
    }
  }
  
  if(any(eigenval<=0)){
    return false;
  }else{
    return true;
  }
  
}



// [[Rcpp::export]]
List initialGen(colvec& y, mat& X, mat& A, int t){
  
  //int n = X.n_rows, p = X.n_cols;
  
  //initial value for betadraw
  colvec betadraw = solve(X,y);
  
  //colvec res = y-X*betadraw;
  
  double sig2draw_y = var(y);
  
  mat Swdraw, Rwdraw = zeros<mat>(t,t) ;
  
  double sigw = 1;
  colvec sig2wdraw(t);
  sig2wdraw.fill(sigw);

  Swdraw = diagmat(sqrt(sig2wdraw));
  Rwdraw.diag().ones();
  
  mat SigWdraw = Swdraw*Rwdraw*Swdraw;
  
 // int dimsigw = SigWdraw.n_rows, 
  int nw = A.n_rows,
      l = nw/t; // # of teachers
  
  colvec muw0 = zeros<colvec>(nw);
  
  mat DiagM = zeros<mat>(l,l);
  DiagM.diag().ones();
  
  mat OmegaW = kron(DiagM, SigWdraw);
  
  colvec wdraw = vectorise(mvrnorm(1,muw0,OmegaW));
  
  colvec gammadraw = solve(A,wdraw);
  
  
  return List::create(Named("betadraw")=betadraw,
                      Named("sig2draw_y")=sig2draw_y,
                      Named("wdraw")=wdraw,
                      Named("sig2wdraw")=sig2wdraw,
                      Named("Swdraw")=Swdraw,
                      Named("Rwdraw")=Rwdraw,
                      Named("SigWdraw")=SigWdraw,
                      Named("gammadraw")=gammadraw);
  
}


// [[Rcpp::export]]
colvec betaUpdate(colvec& y, mat& X, mat& Z, colvec& wdraw, 
                  double sig2draw_y, colvec& B0, mat& sigB){
  
  mat VB;
  colvec muB;
  colvec betadraw;
  
  VB = inv(inv(sigB)+(trans(X)*X)/sig2draw_y);  
  muB = VB*(inv(sigB)*B0+(trans(X)*(y-Z*wdraw))/sig2draw_y);
  
  betadraw = vectorise(mvrnorm(1, muB,VB));
  
  return betadraw;
}
                    



// [[Rcpp::export]]
double sig2yUpdate(colvec& y, mat& X, mat& Z, colvec& wdraw,
                   colvec& betadraw, double nu0_y, double tau0_y){
  
  int n = y.n_elem;
  
  double nu1, tau1, sig2draw_y;
  
  colvec res;
  
  nu1 = (n+2*nu0_y)/2;
  
  res = y-(X*betadraw+Z*wdraw);
  
  tau1 = 0.5*as_scalar(trans(res)*res)+tau0_y;
  
  sig2draw_y = 1/(R::rgamma(nu1,1/tau1));
  
  return sig2draw_y;
  
}



// [[Rcpp::export]]

colvec wUpdate(colvec& y, mat& X, mat& Z, mat& A, int l, colvec& betadraw,
               double sig2draw_y, colvec& gammadraw, mat& Swdraw, mat& Rwdraw){
  
 // int nw = A.n_rows;
 //     l = nw/k;
  
  colvec wdraw;
  mat SigWdraw = Swdraw*Rwdraw*Swdraw;
  
  mat DiagM = zeros<mat>(l,l);
  DiagM.diag().ones();
  
  mat OmegaW = kron(DiagM, SigWdraw);
  
  colvec muW;
  mat VW;
  
  VW = inv((trans(Z)*Z)/sig2draw_y+inv(OmegaW));
  muW = VW*(inv(OmegaW)*(A*gammadraw)+(trans(Z)*(y-X*betadraw))/sig2draw_y);
  
  wdraw = vectorise(mvrnorm(1, muW,VW));
  
  return wdraw;
  
}


// [[Rcpp::export]]

colvec gammaUpdate(colvec& wdraw, mat& A, int l, mat& Swdraw, 
                   mat& Rwdraw, colvec& G0, mat& sigG){
  
  //int nw = A.n_rows;
  //   k = Swdraw.n_rows,l = nw/k;
  
  mat SigWdraw = Swdraw*Rwdraw*Swdraw;
  
  mat DiagM = zeros<mat>(l,l);
  DiagM.diag().ones();
  
  mat OmegaW = kron(DiagM, SigWdraw);
  
  
  mat VG;
  colvec muG;
  colvec gammadraw;
  
  VG = inv(inv(sigG)+(trans(A)*inv(OmegaW)*A));  
  muG = VG*(inv(sigG)*G0+trans(A)*inv(OmegaW)*wdraw);
  
  gammadraw = vectorise(mvrnorm(1, muG,VG));
  
  return gammadraw;
  
}





//[[Rcpp::export]]

double Swpost(colvec& wdraw, mat& A, int t, int l, colvec& gammadraw, 
              mat Swdraw, mat Rwdraw, colvec s0, mat sigS0){
  // sig2wdraw is the variances of w_i 
  // SWdraw is the standard deviation matrix
  //int k = SWdraw.n_rows, 
  
  int nw = wdraw.n_elem;
   
  //mat Swdraw = diagmat(sqrt(sig2wdraw)),
  mat invR = inv(Rwdraw),
      invS = inv(Swdraw),
      invSig = invS*invR*invS;
  
  colvec res;
  
  
  double empSS = 0, logSw, logprior; 
  
  
  for (int j=0; j<nw; j+=t){
    
    colvec res = wdraw.rows(j,j+t-1)-A.rows(j,j+t-1)*gammadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  logprior = as_scalar(-0.5*trans(log(Swdraw.diag())-s0)*inv(sigS0)*(log(Swdraw.diag())-s0));
  
  logSw = as_scalar(-l*log(det(Swdraw))-0.5*empSS)+logprior;
  
  
  return logSw;
  
}


// [[Rcpp::export]]
List SwUpdate(colvec& wdraw, mat& A, int t, int l, colvec& gammadraw, 
              colvec sig2wdraw, mat Swdraw, mat Rwdraw, colvec s0, 
              mat sigS0, colvec sigw_sd, colvec Swcount){
  
  // sig2wdraw  is variances of W_ik
  
  //int k = Swdraw.n_rows;
  colvec sig2wnew = sig2wdraw;
  
  mat Swnew = Swdraw;
  
  double logSwnewpost, logSwdrawpost, accprob;
  
  double logu = log(R::runif(0,1));
  
  
  for(int j = 0;j<t;j++){
    
    sig2wnew = sig2wdraw;
    
    sig2wnew(j) = sig2wdraw(j)+R::rnorm(0,sigw_sd(j));
    
    if(sig2wnew(j)>0){
      Swnew = diagmat(sqrt(sig2wnew));
      logSwnewpost = Swpost(wdraw,A,t,l,gammadraw,Swnew,Rwdraw,s0,sigS0);
    }else{
      logSwnewpost = R_NegInf;
    }
    
    logSwdrawpost = Swpost(wdraw,A,t,l,gammadraw,Swdraw,Rwdraw,s0,sigS0);
    
    accprob = logSwnewpost - logSwdrawpost;
    
    
    if(logu<accprob){
      
      sig2wdraw(j) = sig2wnew(j);
      Swdraw = Swnew;
      Swcount(j) += 1;    
    }     
  }  
 

return List::create(Named("sig2wdraw")=sig2wdraw,
                    Named("Swdraw")=Swdraw,
                    Named("Swcount")=Swcount);

}





//[[Rcpp::export]]

double Rwpost(colvec& wdraw, mat& A, int t, int l, 
              colvec& gammadraw, mat Swdraw, mat Rwdraw){
  
    int  nw = wdraw.n_elem;
      //  k = Swdraw.n_rows, l = nw/k;
      
    mat invR = inv(Rwdraw),
        invS = inv(Swdraw),
        invSig = invS*invR*invS;
      
    colvec res;
      
    double empSS = 0, logRw; 
      
      
    for (int j=0; j<nw; j+=t){
        
      colvec res = wdraw.rows(j,j+t-1)-A.rows(j,j+t-1)*gammadraw;
      empSS += as_scalar(trans(res)*invSig*res);
        
    }
      
  
  logRw = -0.5*l*log(det(Rwdraw))-0.5*empSS;
  
  return logRw;
  
}



// [[Rcpp::export]]
List RwUpdate(colvec& wdraw, mat& A, int t, int l, colvec& gammadraw, 
              mat Swdraw, mat Rwdraw, double division, colvec Rwcount){
   //double rhowdivision,
    int ki, kj;  //k = Rwdraw.n_rows,
    
    mat Rwnew = Rwdraw;
    
    double rhownew, logRwnewpost, logRwdrawpost, accprob;
  
  
    double logu = log(R::runif(0,1));
  
   
   for(ki =0;ki<t-1;ki++){
     
     for(kj=ki+1;kj<t;kj++){
       
       rhownew = Rwdraw(ki,kj)+R::rnorm(0, division);
       
       Rwnew(ki,kj) = rhownew;
       Rwnew(kj,ki) = rhownew;
     }
   }
   
   Rwnew.diag().ones();
   
   if(ispd(Rwnew)){
     logRwnewpost = Rwpost(wdraw,A,t,l,gammadraw,Swdraw,Rwnew);
   }else{
     logRwnewpost = R_NegInf;
   }
    
    logRwdrawpost = Rwpost(wdraw,A,t,l,gammadraw,Swdraw,Rwdraw);
    
    accprob = logRwnewpost - logRwdrawpost;
    
    
    if(logu<accprob){
      
      Rwdraw = Rwnew;
      Rwcount = Rwcount + 1;    
    }     
  
    colvec rhowdraw = corvec(Rwdraw);
  
  
  return List::create(Named("rhowdraw")=rhowdraw,
                      Named("Rwdraw")=Rwdraw,
                      Named("Rwcount")=Rwcount);
  
}




//[[Rcpp::export]]

double logLik(colvec& y, mat& X, mat& Z, mat& A, int t, int l, colvec& betadraw, 
              colvec& wdraw, double sig2draw_y, colvec& gammadraw, mat& Swdraw,
              colvec& sig2wdraw, mat& Rwdraw, colvec& B0, mat& sigB, double nu0_y, 
              double tau0_y, colvec& G0, mat& sigG, colvec s0, mat sigS0){
  
  // k = 2 if w_i = c(w_i1,w_i2), k = 1, if w_i is widiff 
  int n = X.n_rows, p = X.n_cols, q = A.n_cols,
      nw = A.n_rows;
  
  mat SigWdraw = Swdraw*Rwdraw*Swdraw;
  
  mat DiagM = zeros<mat>(l,l);
  DiagM.diag().ones();
  
  mat OmegaW = kron(DiagM, SigWdraw);
  
  
  
  double logf, logy, logbeta,logsig2y, logw, logsigw, loggamma;
  
  
  colvec resy, resw;
  
  resy = y-(X*betadraw+Z*wdraw);
  resw = wdraw-A*gammadraw;
  
  logy = as_scalar(-0.5*n*log(2*M_PI*sig2draw_y)-0.5*trans(resy)*resy/sig2draw_y);
  
  logbeta = as_scalar(-0.5*p*log(2*M_PI)-0.5*log(det(sigB))-0.5*trans(betadraw-B0)*inv(sigB)*(betadraw-B0));
  
  logsig2y = as_scalar(nu0_y*log(tau0_y)-log(Rf_gammafn(nu0_y))-(nu0_y+1)*log(sig2draw_y)-tau0_y/sig2draw_y);
  
  logw = as_scalar(-0.5*nw*log(2*M_PI)-0.5*nw*log(det(SigWdraw))-0.5*trans(resw)*inv(OmegaW)*(resw));
  
  loggamma = as_scalar(-0.5*q*log(2*M_PI)-0.5*log(det(sigG))-0.5*trans(gammadraw-G0)*inv(sigG)*(gammadraw-G0));
  
  logsigw = as_scalar(-0.5*t*log(2*M_PI)-0.5*log(det(sigS0))-0.5*(trans(log(Swdraw.diag())-s0))*inv(sigS0)*(log(Swdraw.diag())-s0));
  
  
  logf = logy+logbeta+logsig2y+logw+loggamma+logsigw;
                        
  return logf;
}





//[[Rcpp::export]]
colvec ypredUpdate(mat& X, mat& Z, colvec& betadraw, colvec& wdraw, double sig2draw_y){
  
  int n = X.n_rows;
  colvec muy = X*betadraw+Z*wdraw;
  
  colvec ypred(n);
  //colvec sigy(n);
  //sigy.fill(sig2draw_y);
  
  //mat Sigy = diagmat(sigy);
  for(int i=0; i<n; i++){
    
    ypred.row(i) = R::rnorm(as_scalar(muy.row(i)),sig2draw_y);
  }
  //colvec ypred = vectorise(mvrnorm(1,muy,Sigy));
                     
  return ypred;                   
                     
}




                    

//[[Rcpp::export]]
List eduSRSmcmcUpdate(colvec& y, mat& X, mat& Z, mat& A, int t, 
                      colvec& B0, mat& sigB, double nu0_y, double tau0_y, 
                      colvec s0, mat sigS0, colvec& G0, mat& sigG, 
                      colvec& sigw_sd, double division, int niter){  //double rhowdivision, 
       //n = X.n_rows, 
       int p = X.n_cols, q = A.n_cols,
           nw = A.n_rows, l = nw/t; 
      // l is # of teachers in both cohorts
      // kt determine the dim of W_i (if (pre-post) k=2, if (wdiff) k = 1) 
       
       colvec betadraw, wdraw, gammadraw, sig2wdraw, rhowdraw;
       //colvec ypreddraw; 
       mat Swdraw, Rwdraw, SigWdraw;                     
       double sig2draw_y;  
       
       //colvec sigWdrawvec; 
       List initialpars, Swout, Rwout;
       
       colvec Swcount = zeros<colvec>(t),
              Rwcount = zeros<colvec>(1);
       double deviance;
       
      // get initial values
      initialpars = initialGen(y,X,A,t);
    
      betadraw = as<colvec>(initialpars["betadraw"]);
      sig2draw_y = as<double>(initialpars["sig2draw_y"]);
      wdraw = as<colvec>(initialpars["wdraw"]);
      sig2wdraw = as<colvec>(initialpars["sig2wdraw"]);
      Swdraw = as<mat>(initialpars["Swdraw"]);
      Rwdraw = as<mat>(initialpars["Rwdraw"]);
      gammadraw = as<colvec>(initialpars["gammadraw"]);
       
      rhowdraw = corvec(Rwdraw);
      int nrho = rhowdraw.n_elem;
      
     
      mat betamat = zeros<mat>(niter,p),
          Wmat = zeros<mat>(niter,nw),
          gammamat = zeros<mat>(niter,q),
          sig2wmat = zeros<mat>(niter,t),
          rhowmat = zeros<mat>(niter,nrho);
      
        //  ypredmat = zeros<mat>(niter,n);
      colvec sig2ymat = zeros<colvec>(niter),
             devmat = zeros<colvec>(niter);
      
      
      
      for(int iter=0; iter<niter; iter++){
        
        betadraw = betaUpdate(y,X,Z,wdraw,sig2draw_y,B0,sigB);
        betamat.row(iter) = trans(betadraw);
        
        sig2draw_y = sig2yUpdate(y,X,Z,wdraw,betadraw,nu0_y,tau0_y);
        sig2ymat(iter) = sig2draw_y;
        
        wdraw = wUpdate(y,X,Z,A,l,betadraw,sig2draw_y, 
                        gammadraw,Swdraw,Rwdraw);
        Wmat.row(iter) = trans(wdraw);
        
        gammadraw = gammaUpdate(wdraw,A,l,Swdraw,Rwdraw,G0,sigG); 
        gammamat.row(iter) = trans(gammadraw);
                   
        Swout = SwUpdate(wdraw,A,t,l,gammadraw,sig2wdraw,Swdraw,
                         Rwdraw,s0,sigS0,sigw_sd,Swcount);
        sig2wdraw = as<colvec>(Swout["sig2wdraw"]);
        Swdraw = as<mat>(Swout["Swdraw"]);
        Swcount = as<colvec>(Swout["Swcount"]);
        sig2wmat.row(iter) = trans(sig2wdraw);
         
        
        Rwout = RwUpdate(wdraw,A,t,l,gammadraw,Swdraw,Rwdraw,division,Rwcount);
        rhowdraw = as<colvec>(Rwout["rhowdraw"]);
        Rwdraw = as<mat>(Rwout["Rwdraw"]);
        Rwcount = as<colvec>(Rwout["Rwcount"]);
        rhowmat.row(iter) = rhowdraw;
        
        
        deviance = -2*logLik(y, X, Z, A, t, l, betadraw, wdraw, sig2draw_y, 
                             gammadraw, Swdraw, sig2wdraw, Rwdraw, B0, sigB, 
                             nu0_y, tau0_y, G0, sigG, s0, sigS0);
        devmat(iter) = deviance;
        
   
      }
       
       
      return List::create(Named("betamat")=betamat,
                          Named("sig2ymat")=sig2ymat,
                          Named("Wmat")=Wmat,
                          Named("gammamat")=gammamat,
                          Named("sig2wmat")=sig2wmat,
                          Named("rhowmat")=rhowmat,
                          Named("devmat")=devmat,
                          Named("Swcount")=Swcount,
                          Named("Rwcount")=Rwcount); 
}


//================================================================================================


//[[Rcpp::export]]
List eduadaptiveUpdate(colvec& y, mat& X, mat& Z, mat& A, int t, colvec& B0, mat& sigB, 
                       double nu0_y, double tau0_y, colvec s0, mat sigS0, colvec& G0,
                       mat& sigG, colvec& sigw_sd, double division, int nadap, int niter){  //double rhowdivision, 
  
  
  const double lowaccrate = 0.2, upaccrate = 0.5;
  
  //lowaccrate = 0.234; 
  //upaccrate = 0.44;
  
  
  double lowacc_count = lowaccrate*nadap, 
    upacc_count = upaccrate*nadap;
  
  //n = X.n_rows, 
  int p = X.n_cols, q = A.n_cols,
    nw = A.n_rows, l = nw/t; 
  // l is # of teachers in both cohorts
  // kt determine the dim of W_i (if (pre-post) k=2, if (wdiff) k = 1) 
  
  colvec betadraw, wdraw, gammadraw, sig2wdraw, rhowdraw;
  //colvec ypreddraw; 
  mat Swdraw, Rwdraw, SigWdraw;                     
  double sig2draw_y;  
  
  //colvec sigWdrawvec; 
  List initialpars, Swout, Rwout;
  
  colvec Swcount = zeros<colvec>(t),
         Rwcount = zeros<colvec>(1);
  
  NumericVector count = wrap(join_cols(Swcount,Rwcount));
  
  
  
  double deviance;
  
  // get initial values
  initialpars = initialGen(y,X,A,t);
  
  betadraw = as<colvec>(initialpars["betadraw"]);
  sig2draw_y = as<double>(initialpars["sig2draw_y"]);
  wdraw = as<colvec>(initialpars["wdraw"]);
  sig2wdraw = as<colvec>(initialpars["sig2wdraw"]);
  Swdraw = as<mat>(initialpars["Swdraw"]);
  Rwdraw = as<mat>(initialpars["Rwdraw"]);
  gammadraw = as<colvec>(initialpars["gammadraw"]);
  
  rhowdraw = corvec(Rwdraw);
  int nrho = rhowdraw.n_elem;
  
  
  mat betamat = zeros<mat>(niter,p),
    Wmat = zeros<mat>(niter,nw),
    gammamat = zeros<mat>(niter,q),
    sig2wmat = zeros<mat>(niter,t),
    rhowmat = zeros<mat>(niter,nrho);
  
  //  ypredmat = zeros<mat>(niter,n);
  colvec sig2ymat = zeros<colvec>(niter),
    devmat = zeros<colvec>(niter);
  
  
  
  double testno = 0;
  
  while(is_true(any(count<lowacc_count)) || is_true(any(count>upacc_count))){
    
    testno += 1;
    cout<<"test "<< testno<<endl;
    
    for(int iadap=0; iadap<nadap; iadap++){
      
      betadraw = betaUpdate(y,X,Z,wdraw,sig2draw_y,B0,sigB);
     
      sig2draw_y = sig2yUpdate(y,X,Z,wdraw,betadraw,nu0_y,tau0_y);
     
      wdraw = wUpdate(y,X,Z,A,l,betadraw,sig2draw_y, 
                      gammadraw,Swdraw,Rwdraw);
      
      gammadraw = gammaUpdate(wdraw,A,l,Swdraw,Rwdraw,G0,sigG); 
      
      Swout = SwUpdate(wdraw,A,t,l,gammadraw,sig2wdraw,Swdraw,
                       Rwdraw,s0,sigS0,sigw_sd,Swcount);
      sig2wdraw = as<colvec>(Swout["sig2wdraw"]);
      Swdraw = as<mat>(Swout["Swdraw"]);
      Swcount = as<colvec>(Swout["Swcount"]);
   
      
      Rwout = RwUpdate(wdraw,A,t,l,gammadraw,Swdraw,Rwdraw,division,Rwcount);
      rhowdraw = as<colvec>(Rwout["rhowdraw"]);
      Rwdraw = as<mat>(Rwout["Rwdraw"]);
      Rwcount = as<colvec>(Rwout["Rwcount"]);
    } 
      
      count = wrap(join_cols(Swcount,Rwcount));
      
      for(int i=0; i<t; i++){
        
        if(Swcount(i)<lowacc_count && sigw_sd(i)>(0.05*sigw_sd(i))){
          sigw_sd(i) -= 0.05*sigw_sd(i);
        }
        
        if(Swcount(i)>upacc_count){
          sigw_sd(i) += 0.05*sigw_sd(i);
        }
        
      }
      
      if(as_scalar(Rwcount)<lowacc_count && as_scalar(division)>as_scalar(0.1*division)){
        division -= 0.1*division;
      }
      
      if(as_scalar(Rwcount)>upacc_count){
        division += 0.1*division;
      }
      
      
      cout<<Swcount<<endl;
      cout<<sigw_sd<<endl;
      cout<<Rwcount<<endl;
      cout<<division<<endl;
      
      Swcount.zeros();
      Rwcount.zeros();
    
  }
  
  cout<<"pass"<<endl;
      
  
  
  for(int iter=0; iter<niter; iter++){
    
    betadraw = betaUpdate(y,X,Z,wdraw,sig2draw_y,B0,sigB);
    betamat.row(iter) = trans(betadraw);
    
    sig2draw_y = sig2yUpdate(y,X,Z,wdraw,betadraw,nu0_y,tau0_y);
    sig2ymat(iter) = sig2draw_y;
    
    wdraw = wUpdate(y,X,Z,A,l,betadraw,sig2draw_y, 
                    gammadraw,Swdraw,Rwdraw);
    Wmat.row(iter) = trans(wdraw);
    
    gammadraw = gammaUpdate(wdraw,A,l,Swdraw,Rwdraw,G0,sigG); 
    gammamat.row(iter) = trans(gammadraw);
    
    Swout = SwUpdate(wdraw,A,t,l,gammadraw,sig2wdraw,Swdraw,
                     Rwdraw,s0,sigS0,sigw_sd,Swcount);
    sig2wdraw = as<colvec>(Swout["sig2wdraw"]);
    Swdraw = as<mat>(Swout["Swdraw"]);
    Swcount = as<colvec>(Swout["Swcount"]);
    sig2wmat.row(iter) = trans(sig2wdraw);
    
    
    Rwout = RwUpdate(wdraw,A,t,l,gammadraw,Swdraw,Rwdraw,division,Rwcount);
    rhowdraw = as<colvec>(Rwout["rhowdraw"]);
    Rwdraw = as<mat>(Rwout["Rwdraw"]);
    Rwcount = as<colvec>(Rwout["Rwcount"]);
    rhowmat.row(iter) = rhowdraw;
    
    
    deviance = -2*logLik(y, X, Z, A, t, l, betadraw, wdraw, sig2draw_y, 
                         gammadraw, Swdraw, sig2wdraw, Rwdraw, B0, sigB, 
                         nu0_y, tau0_y, G0, sigG, s0, sigS0);
    devmat(iter) = deviance;
    
    
  }
  
  
  return List::create(Named("betamat")=betamat,
                      Named("sig2ymat")=sig2ymat,
                      Named("Wmat")=Wmat,
                      Named("gammamat")=gammamat,
                      Named("sig2wmat")=sig2wmat,
                      Named("rhowmat")=rhowmat,
                      Named("devmat")=devmat,
                      Named("Swcount")=Swcount,
                      Named("Rwcount")=Rwcount); 
}


