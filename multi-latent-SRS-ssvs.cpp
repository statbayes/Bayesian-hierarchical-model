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


//=========================================================================
// get upper/lower triangular without diagonal (get correlations vector)
//=========================================================================

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


//=====================================================================
// generate from multivariate normal
//=====================================================================

// [[Rcpp::export]]
mat mvrnorm(int n, colvec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return (repmat(mu, 1, n).t() + Y * chol(sigma));
}


const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
colvec dmvnorm(mat x, rowvec mu,  mat sigma, 
               bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  colvec out(n);
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  colvec z;
  
  for (int i=0; i < n; i++) {
    z = rooti * trans( x.row(i) - mu) ;    
    out(i)  = constants - 0.5 * sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}


//===========================================================================
// generate from wishart and inverse wishart
//===========================================================================

// wishart and inverse wishart

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



//=================================================================
// Generate initial values
//=================================================================
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



//===================================================================
// Update parameters of fixed effects for students 
// (intercept,lowinc, eth, lowinc*eth)
//===================================================================


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


//====================================================================
// Update variance for students' scores
//====================================================================

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

//===================================================================
// Update pre- and post- teachers' latent effects
//===================================================================

// [[Rcpp::export]]

colvec wUpdate(colvec& y, mat& X, mat& Z, mat& A, int l, colvec& betadraw,
               double sig2draw_y, colvec& gammadraw, mat& Swdraw, mat& Rwdraw){
  
  //int dimsigw = SigWdraw.n_rows, 
  //  nrowA = A.n_rows,
  //  nw = nrowA/dimsigw;
  
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


//=============================================
// Update Sw
//===============================================


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


//=================================================
// Update Rwdraw
//=================================================

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












//===============================================================
// Update parameters of teachers' instruments 
//(predictors for teachers' latent effects)
//===============================================================

// [[Rcpp::export]]

colvec gammaUpdate(colvec& wdraw, mat& A, int l, mat& Swdraw, 
                   mat& Rwdraw, colvec& G0, mat& sigG){
  
  
  //mat sigG = SigGUpdate(sigmag0, delta, c2, t2);
  
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


//=================================================================================
// Update prior var-covariance matrix sigG, 
// delta_ij*N(0,c^2*t^2)+(1-delta_ij)*N(0,t^2)
// if the variable is always selected, set prior prob of gamma_ij,  p_ij=1
//=================================================================================

//[[Rcpp::export]]

mat priorSigUpdate(colvec& selectindex, colvec& cidx, colvec& tidx, mat& R0){
  // c and t are the scale (c^2) and precision(tidx^2)
  
  int p = selectindex.n_elem;
  
  colvec a = ones<colvec>(p);
  
  mat V;
  
  for(int j=0; j<p; j++){
    
    if(selectindex(j)==1){
      a(j) = cidx(j);
    }
  }
  
  colvec dvec = a%tidx;
  
  mat D_idx = diagmat(dvec);
  
  V = D_idx*R0*D_idx;
  
  return V;
}

//==========================================================================
// likelihood of gammadraw, ~ N_p(0_p, sigG)
// sigG = priorSigUpdate(deltadraw,cidx_g,tidx_g,RG0),
//==========================================================================

//[[Rcpp::export]]

colvec repvec(colvec& a, double m){
  
  NumericVector aa = as<NumericVector>(wrap(a)),
    v = rep(aa,m);
  return as<colvec>(v);
}




//[[Rcpp::export]]

double gammalikelihood(colvec& gammadraw, colvec& deltadraw, colvec& c_g, 
                       colvec& t_g, mat& RG0, mat& IG){
  
  
  mat HG = priorSigUpdate(deltadraw,c_g,t_g,RG0); 
  
  mat SigG = kron(IG, HG);
  
  //cout<<SigA<<endl;
  
  int pk = SigG.n_cols;
  
  mat invS = inv(SigG);
  
  double loggamma = as_scalar(-0.5*pk*log(2*M_PI)-0.5*log(det(SigG))-0.5*trans(gammadraw)*invS*gammadraw);
  
  return exp(loggamma);
  
} 


//[[Rcpp::export]]

double gammalikelihood1(colvec& gammadraw, colvec& deltadraw, 
                        colvec& c_g, colvec& t_g, mat& RG0){
  
  
  mat SigG = priorSigUpdate(deltadraw,c_g,t_g,RG0); 
  
  //mat SigG = kron(IG, HG);
  
  //cout<<SigA<<endl;
  
  int pk = SigG.n_cols;
  
  mat invS = inv(SigG);
  
  double loggamma = as_scalar(-0.5*pk*log(2*M_PI)-0.5*log(det(SigG))-0.5*trans(gammadraw)*invS*gammadraw);
  
  return exp(loggamma);
  
} 


//======================================================================
// Update selection index , deltadraw
//======================================================================

//[[Rcpp::export]]

colvec deltaUpdate(colvec& gammadraw, colvec& deltadraw, colvec& deltaprob0, 
                   colvec& c_g, colvec& t_g, mat& RG0, mat IG){
  
  // parsdraw = gammadraw
  // prob0 = prior probabilities of selectindex(deltadraw)
  // set prob0[j] = 1 if we always include the variable,
  // e.g. always include 'intercept' then prob0 = c(1,0.5,0.5,1,0.5,...)
  
  
  int np = deltadraw.n_elem; // length(parsdraw) = length(prob0)
  
  double p0, p1, postp_j; //p0 with indexdraw[j]=0, p1 with indexdraw[1]
  
  
  for(int j=0; j<np; j++){
    
    deltadraw(j) = 0;
    p0 = (1-deltaprob0(j))*gammalikelihood(gammadraw,deltadraw,c_g,t_g,RG0,IG);
    
    deltadraw(j) = 1;   
    p1 = deltaprob0(j)*gammalikelihood(gammadraw,deltadraw,c_g,t_g,RG0,IG);
    
    postp_j = p1/(p0+p1);
    
    deltadraw(j) = R::rbinom(1,postp_j);
    
  }
  
  return deltadraw;
}



//[[Rcpp::export]]

colvec deltaUpdate1(colvec& gammadraw, colvec& deltadraw, colvec& deltaprob0, 
                    colvec& c_g, colvec& t_g, mat& RG0){
  
  // parsdraw = gammadraw
  // prob0 = prior probabilities of selectindex(deltadraw)
  // set prob0[j] = 1 if we always include the variable,
  // e.g. always include 'intercept' then prob0 = c(1,0.5,0.5,1,0.5,...)
  
  
  int np = deltadraw.n_elem; // length(parsdraw) = length(prob0)
  
  double p0, p1, postp_j; //p0 with indexdraw[j]=0, p1 with indexdraw[1]
  
  
  for(int j=0; j<np; j++){
    
    deltadraw(j) = 0;
    p0 = (1-deltaprob0(j))*gammalikelihood1(gammadraw,deltadraw,c_g,t_g,RG0);
    
    deltadraw(j) = 1;   
    p1 = deltaprob0(j)*gammalikelihood1(gammadraw,deltadraw,c_g,t_g,RG0);
    
    postp_j = p1/(p0+p1);
    
    deltadraw(j) = R::rbinom(1,postp_j);
    
  }
  
  return deltadraw;
}





//===========================================================================
// MCMC update with Variable selection (SSVS)
//===========================================================================


//[[Rcpp::export]]
List SRSssvsUpdate(colvec& y, mat& X, mat& Z, mat& A, int t, colvec& B0, mat& sigB, 
                double nu0_y, double tau0_y, colvec s0, mat sigS0, colvec& G0, 
                colvec& deltaprob0, colvec& c_g, colvec& t_g, mat& RG0, 
                colvec& sigw_sd, double division, int niter){
  
  //colvec& sigmag0, 
  //int n=X.n_rows;
  int p = X.n_cols, q = A.n_cols, 
    nw = A.n_rows, l = nw/t; // k is the time, k=2 pre-post
  int qd = q/t;
  
  colvec betadraw, wdraw, gammadraw, sig2wdraw, rhowdraw;
  //colvec ypreddraw; 
  mat Swdraw, Rwdraw, SigWdraw;                
  double sig2draw_y;  
  //double deviance;
  
  //colvec deltadraw = ones<colvec>(q-k);
  colvec deltadraw = ones<colvec>(qd);
 
  colvec Swcount = zeros<colvec>(t),
         Rwcount = zeros<colvec>(1);
 
  List initialpars, Swout, Rwout;
  
  // for(int j=0;j<q-1;j++){
  //   deltadraw(j) = R::rbinom(1,p0);
  // }
  
  
  
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
  

  //mat sigG = SigGUpdate(sigmag0, deltadraw, c2, t2);
  mat HG = priorSigUpdate(deltadraw,c_g,t_g,RG0),
      IG(t,t);
  IG.diag().ones();
  mat sigG = kron(IG,HG);
  
  
  mat betamat = zeros<mat>(niter,p),
    Wmat = zeros<mat>(niter,nw),
    sig2wmat = zeros<mat>(niter,t),
    rhowmat = zeros<mat>(niter,nrho),
    gammamat = zeros<mat>(niter,q),
    deltamat = zeros<mat>(niter,qd);
  // deltamat = zeros<mat>(niter,q-k)
  //mat ypredmat = zeros<mat>(niter,n);
  colvec sig2ymat = zeros<colvec>(niter);// devmat(niter);
  
  
  
  for(int iter=0; iter<niter; iter++){
    
    
    betadraw = betaUpdate(y, X, Z, wdraw, sig2draw_y, B0, sigB);
    betamat.row(iter) = trans(betadraw);
    
    sig2draw_y = sig2yUpdate(y,X,Z,wdraw,betadraw,nu0_y,tau0_y);
    sig2ymat[iter] = sig2draw_y;
    
    wdraw = wUpdate(y,X,Z,A,l,betadraw,sig2draw_y,gammadraw,Swdraw,Rwdraw); 
    Wmat.row(iter) = trans(wdraw);
    
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
    
    
    gammadraw = gammaUpdate(wdraw,A,l,Swdraw,Rwdraw,G0,sigG);
    gammamat.row(iter) = trans(gammadraw);
    
    deltadraw = deltaUpdate(gammadraw, deltadraw, deltaprob0, c_g, t_g, RG0, IG);
    deltamat.row(iter) = trans(deltadraw);
    
    HG = priorSigUpdate(deltadraw, c_g, t_g, RG0);
    sigG = kron(IG,HG);
    
    
    // deviance = -2*logLik(y, X, Z, A, betadraw, wdraw, sig2draw_y, 
    //                   gammadraw, SigWdraw, B0, sigB, nu0_y, tau0_y,
    //                   nuw, Sw, G0, sigG);
    //ypreddraw = ypredUpdate(X, Z, betadraw, wdraw, sig2draw_y);
    
    //ypredmat.row(iter) = trans(ypreddraw);
    // devmat[iter] = deviance;
  }
  
  
  return List::create(Named("betamat")=betamat,
                      Named("sig2ymat")=sig2ymat,
                      Named("Wmat")=Wmat,
                      Named("sig2wmat")=sig2wmat,
                      Named("rhowmat")=rhowmat,
                      Named("gammamat")=gammamat,
                      Named("deltamat")=deltamat,
                      Named("Swcount")=Swcount,
                      Named("Rwcount")=Rwcount);
  //Named("devmat")=devmat); 
}



//[[Rcpp::export]]
List SSVSUpdate1(colvec& y, mat& X, mat& Z, mat& A, int t, colvec& B0, mat& sigB, 
                 double nu0_y, double tau0_y, colvec s0, mat sigS0, colvec& G0, 
                 colvec& deltaprob0, colvec& c_g, colvec& t_g, mat& RG0, 
                 colvec& sigw_sd, double division, int niter){
  
  //colvec& sigmag0, 
  //int n=X.n_rows;
  int p = X.n_cols, q = A.n_cols, 
    nw = A.n_rows, l = nw/t; // k is the time, k=2 pre-post
  int qd = q/t;
  
  colvec betadraw, wdraw, gammadraw, sig2wdraw, rhowdraw;
  //colvec ypreddraw; 
  mat Swdraw, Rwdraw, SigWdraw;                
  double sig2draw_y;  
  //double deviance;
  
  //colvec deltadraw = ones<colvec>(q-k);
  colvec deltadraw = ones<colvec>(qd);
  
  colvec Swcount = zeros<colvec>(t),
    Rwcount = zeros<colvec>(1);
  
  List initialpars, Swout, Rwout;
  
  // for(int j=0;j<q-1;j++){
  //   deltadraw(j) = R::rbinom(1,p0);
  // }
  
  
  
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
  
  
  //mat sigG = SigGUpdate(sigmag0, deltadraw, c2, t2);
  mat HG = priorSigUpdate(deltadraw,c_g,t_g,RG0),
    IG(t,t);
  IG.diag().ones();
  mat sigG = kron(IG,HG);
  
  
  mat betamat = zeros<mat>(niter,p),
    Wmat = zeros<mat>(niter,nw),
    sig2wmat = zeros<mat>(niter,t),
    rhowmat = zeros<mat>(niter,nrho),
    gammamat = zeros<mat>(niter,q),
    deltamat = zeros<mat>(niter,qd);
  // deltamat = zeros<mat>(niter,q-k)
  //mat ypredmat = zeros<mat>(niter,n);
  colvec sig2ymat = zeros<colvec>(niter);// devmat(niter);
  
  
  for(int iter=0; iter<niter; iter++){
    
    
    betadraw = betaUpdate(y, X, Z, wdraw, sig2draw_y, B0, sigB);
    betamat.row(iter) = trans(betadraw);
    
    sig2draw_y = sig2yUpdate(y,X, Z, betadraw,wdraw, nu0_y, tau0_y);
    sig2ymat[iter] = sig2draw_y;
    
    wdraw = wUpdate(y,X,Z,A,l,betadraw,sig2draw_y,gammadraw,Swdraw,Rwdraw); 
    Wmat.row(iter) = trans(wdraw);
    
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
    
    
    gammadraw = gammaUpdate(wdraw,A,l,Swdraw,Rwdraw,G0,sigG);
    gammamat.row(iter) = trans(gammadraw);
    
    deltadraw = deltaUpdate1(gammadraw, deltadraw, deltaprob0, c_g, t_g, RG0);
    deltamat.row(iter) = trans(deltadraw);
    
    sigG = priorSigUpdate(deltadraw, c_g, t_g, RG0);
    //sigG = kron(IG,HG);
    
    
    // deviance = -2*logLik(y, X, Z, A, betadraw, wdraw, sig2draw_y, 
    //                   gammadraw, SigWdraw, B0, sigB, nu0_y, tau0_y,
    //                   nuw, Sw, G0, sigG);
    //ypreddraw = ypredUpdate(X, Z, betadraw, wdraw, sig2draw_y);
    
    //ypredmat.row(iter) = trans(ypreddraw);
    // devmat[iter] = deviance;
  }
  
  
  return List::create(Named("betamat")=betamat,
                      Named("sig2ymat")=sig2ymat,
                      Named("Wmat")=Wmat,
                      Named("sig2wmat")=sig2wmat,
                      Named("rhowmat")=rhowmat,
                      Named("gammamat")=gammamat,
                      Named("deltamat")=deltamat);
  //Named("devmat")=devmat); 
}


