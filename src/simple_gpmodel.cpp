#include "simple_gpmodel.h"

//p1 p2: position of first and second paraeter
//k==0: kernel value; kernel_K
//k==1 partial deriv w.r.t. theta_1, kernel_partialK (output)
//k==2: partial deriv w.r.t theta_2, kernel_partialK (output)
//k(xx, yy)
//If K is symmetric, save half computation time
void simple_gpModel::kernel_SE(int p1, int p2,
                               int k, vec & xx, vec & yy,
                               bool symmetric){
  int i,j, jStart;
  int nK1, nK2;
  bool zero_out;
  nK1 = xx.n_elem;
  nK2 = yy.n_elem;
  double tmpE,tmpE_deriv, tmpK, tmpDiff;
  kernel_K = mat(nK1, nK2, fill::zeros);
  kernel_partialK = kernel_K;
  
  double eTheta1,eTheta2;
  eTheta1 = exp(theta(p1));
  eTheta2 = exp(theta(p2));
  
  jStart = 0;
  for(i=0; i< nK1; i++){
    if(symmetric) jStart = i;
    for(j=jStart; j<nK2; j++){
      tmpDiff = xx(i) - yy(j);
      zero_out = false;
      //if(abs(tmpDiff)>1500) zero_out = true;
      //long term trend
      tmpE = tmpDiff * tmpDiff /(2*eTheta2*eTheta2);
      if(k==0){
        tmpK = eTheta1*eTheta1*exp(-tmpE);
        kernel_K(i,j) = tmpK;
        if(zero_out) kernel_K(i,j) = 0;
        if(symmetric) kernel_K(j,i) = kernel_K(i,j);
      }else if(k==1){
        tmpK = 2*eTheta1*exp(-tmpE)*eTheta1;
        kernel_partialK(i,j) = tmpK;
        if(zero_out) kernel_partialK(i,j) = 0;
        if(symmetric) kernel_partialK(j,i) =
          kernel_partialK(i,j);
      }else if(k==2){
        tmpK = eTheta1*eTheta1*exp(-tmpE);
        tmpE_deriv = 2*tmpE; // /eTheta2*eTheta2
        tmpK *= tmpE_deriv;
        kernel_partialK(i,j) = tmpK;
        if(zero_out) kernel_partialK(i,j) = 0;
        if(symmetric) kernel_partialK(j,i) =
          kernel_partialK(i,j);
      }
    }// for j
  }// for i
}// kernel_SE


//compute the variance/covariance matrix: k(xx, yy)
//k==-1: training data, used in prediction
//k==0: training data, theta estimation
//k==1: cov(new, train) --> gp_K21
//k==2: var(new) --> gp_K22
//
// includeBeta, the sinusoid component, MIRA/SRV, k must be 0
// includeGP, the GP component, pure sinusoid, k must be 0
// k==0 for pure likelihood evaluation
void simple_gpModel::compute_CovM(vec & xx, vec & yy, int k,
                                  bool includeBeta=true,
                                  bool includeGP=true){
  // compute kernel
  bool symmetric = false;
  if(k<=0) symmetric = true;
  
  int nK1,nK2;
  nK1 = xx.n_elem;
  nK2 = yy.n_elem;
  gp_Kc = mat(nK1, nK2, fill::zeros);

  
  
  // for training, include measurement error
  // else, leave early
  if(k<=0){
    for(int i=0; i<nObs; i++){
      gp_Kc(i, i)  += mag_sigma(i)*mag_sigma(i)
      + 0.001 ;  //relaxed
    }
  }else if (k==1){
    kernel_SE(0, 1, 0, xx, yy, symmetric);
    gp_K21 = kernel_K;
    return;
  }else if(k==2){
    kernel_SE(0, 1, 0, xx, yy, symmetric);
    gp_K22 = kernel_K;
    return;
  }
  
 
  //include GP kernel
  if(includeGP){
    kernel_SE(0, 1, 0, xx, yy, symmetric);
    gp_Kc += kernel_K;
  }else{
     double eTheta1 = exp(theta(0));
     eTheta1 *= eTheta1;
    gp_Kc.diag() += eTheta1;
  }


  //if Mira model, add HtBH
  //else only add overall mean uncertainty
  if(includeBeta)
    gp_Ky = gp_Kc + HtBH;
  else
    gp_Ky = gp_Kc + sigmamSq;

  //compute chol, for inverse
  flag_inv = true;
  flag_inv = chol(gp_Ky_chol,gp_Ky);
  if(flag_inv){
    vec gp_Ky_diag = gp_Ky_chol.diag();
    gp_Ky_logdet = sum(log(gp_Ky_diag));
    gp_Ky_logdet = 2 * gp_Ky_logdet;
    //inverse
    alpha = solve(trimatl(gp_Ky_chol.t()), gp_y);
    alpha = solve(trimatu(gp_Ky_chol), alpha);
    if(k==-1) gp_Kc_inv = gp_Kc.i();
  }// flag_inv
  
}// compute_CovM



//copmute minus log likelihood of the data
// this the joint density of y and theta
// includeBeta = false, model for SRV
// includeGP = false, model for O-rich?? or C-rich
double simple_gpModel::gp_mLoglik(vec theta_,
                                  bool includeBeta = true,
                                  bool includeGP = true){
  double mLoglik;
  //gp_setTheta(theta_);
  theta = theta_;
  
  compute_CovM(MJD, MJD, 0, includeBeta, includeGP);
  
  if(flag_inv){
    //&& abs(theta(0))<6  && theta(1)>2 && theta(1)<500
    //log_det(mLoglik, signL, gp_Ky);
    mLoglik =0.5*gp_Ky_logdet +
      0.5*nObs*log(2*PI)+
      0.5* as_scalar(gp_y.t()*alpha);         
  }else{
    mLoglik = largeV;
  }    
  return mLoglik;
}// gp_mLoglik



//compute the partial derivative of the minus log likelihood
vec simple_gpModel::gp_DmLoglik(vec theta_){
  
  vec partialDeri(2,fill::zeros);
  //if(!flag_inv) return partialDeri;
  theta = theta_;
  
  //compute_CovM(MJD, MJD, 0); // save a little time by Optim
  
  mat tempI = eye(nObs,nObs);
  gp_Ky_inv = solve(trimatl(gp_Ky_chol.t()), tempI);
  gp_Ky_inv = solve(trimatu(gp_Ky_chol), gp_Ky_inv);
  
  mat prodTmp;
  prodTmp = alpha*alpha.t() - gp_Ky_inv;
  //kernel
  int j;
  for(j=0; j<2;j++){
    kernel_SE(0,1, j+1, MJD, MJD, true);
    partialDeri(j) = -0.5*trace(prodTmp*kernel_partialK);
  }
  return partialDeri;
}// gp_DmLoglik





//set theta 
void simple_gpModel::set_theta(NumericVector theta_){
  theta = vec(theta_);
}


//The design matrix: n-by-3
//(1,cos(xx), sin(xx))
mat simple_gpModel::gp_get_Ht(vec xx){
  int nK = xx.n_elem;
  mat Ht_tmp = mat(nK, 3, fill::ones);
  xx = xx*2*PI * freq;
  Ht_tmp(span::all, 1) = cos(xx);
  Ht_tmp(span::all, 2) = sin(xx);
  return Ht_tmp;
}// gp_get_Ht

void simple_gpModel::set_freq(double fff, double shift){
  freq = fff;
  //new design matrix after changing period
  Ht = gp_get_Ht(MJD);
  HtBH =  Ht*Sigma_gamma*Ht.t();
}//set_freq


//predict at Xnew (t*)
//Posterior distribution of gamma and Ynew
List simple_gpModel::predict(NumericVector Xnew_){
  Xnew = vec(Xnew_);
  
  vec resid_hat;
  //mat resid_Omega,resid_Sigma;
  //double mLoglik_resid = 0, signL;
  
  //-1 cov of training data, gp_Kc_inv is also computed
  compute_CovM(MJD, MJD, -1);
  
  
  compute_CovM(Xnew, MJD, 1);
  compute_CovM(Xnew, Xnew, 2);
  
  //get the design matrix
  mat Ht_star = gp_get_Ht(Xnew);
  
  //posterior mean of gamma and new Y
  vec gamma_bar, predY;
  //posterior covariance matrix of gamma and new Y
  mat tmp1, tmp2, tmp_r, pcov_predY,pcov_gamma;
  
  tmp1 = Sigma_gamma.i();
  tmp1 += Ht.t()*gp_Kc_inv*Ht;
  tmp2 = Sigma_gamma.i()*gamma0+ Ht.t()*gp_Kc_inv*(gp_y+gp_m);
  
  pcov_gamma = inv(tmp1);
  gamma_bar = pcov_gamma*tmp2;
  
  resid_hat = gp_y+ gp_m- Ht*gamma_bar;
  predY = Ht_star*gamma_bar + 
    gp_K21 * gp_Kc_inv * resid_hat;
  tmp_r = Ht_star.t() - Ht.t() * gp_Kc_inv *gp_K21.t();
  pcov_predY = gp_K22 - gp_K21 * gp_Kc_inv * gp_K21.t();
  pcov_predY += tmp_r.t() * pcov_gamma * tmp_r;
  //    log_det(mLoglik_resid, signL, resid_Sigma);
  //mLoglik_resid =0.5*mLoglik_resid +
  //	0.5* as_scalar(resid_hat.t()*resid_Omega*resid_hat);
  
  return List::create(
    Named("Xnew") = Xnew,
    Named("predy")=predY,
    Named("predy_cov") = pcov_predY,
    Named("gamma_bar")=gamma_bar,
    Named("gamma_cov") = pcov_gamma,
    Named("Ht_star") = Ht_star
  );
}



simple_gpModel::simple_gpModel(NumericVector MJD_,
                               NumericVector mag_, NumericVector error_):
                               varStar(MJD_, mag_, error_)
{
  gp_m = 15.62 + 6.2;
  sigmamSq = 100;
  sigmabSq= 1;
  fmin = 0.0005;
  fmax = 0.01;
  fdelta = 0.00001;
  
  //priors for gamma
  gamma0 = vec(3,fill::zeros);
  gamma0(0) = gp_m;
  Sigma_gamma = mat(3, 3, fill::zeros);
  Sigma_gamma(0,0) = sigmamSq;
  Sigma_gamma(1,1) = sigmabSq;
  Sigma_gamma(2,2) = sigmabSq;
  // remove the overall prior mean
  gp_y = mag - gp_m;
  
  largeV = 1.0e30;
}// constructor





// the optimization method for frequency estimation
mat simple_gpModel::freq_est(){
  double ftrial;
  int fstep;
  
  fstep = floor((fmax-fmin)/fdelta);
  
  
  mat spc = mat(fstep,4,fill::zeros);
  vec theta0(2,fill::zeros),theta0Opt(2,fill::zeros);
  mat H0,H0Opt;
  double yopt,cloglik;
  
  
  ftrial = fmax;
  int i, j,h;
  int k=0, checkL = 90;
  
  h = 0;
  // length(seq(-2,15.9,by=1))* length(seq(-2,2.9,by=1))
  
  mat init = mat(2,checkL,fill::zeros);
  vec initY(checkL,fill::zeros);
  uword minRowI;
  
  for(i=-2;i<3;i=i+1)
    for(j = -2;j<16;j=j+1){
      init(0,h) = (double) i;
      init(1,h) = (double) j;
      h++;
    }
    
    
    while(ftrial>=fmin){
      set_freq(ftrial,0);
      
      cloglik = gp_mLoglik(theta0Opt,true,true);
      
      // evaluate the whole surface
      for(h=0;h<checkL;h++){
        initY(h) = gp_mLoglik(init.col(h),true,true);
      }
      
      //if the current opt point is not sub-optimal
      //then restart
      if(cloglik>initY.min(minRowI)){
        H0 = eye<mat>(2,2);
        theta0 = init.col(minRowI);
      }else{
        theta0 = theta0Opt;
        H0 = H0Opt;
      }
      
      
      yopt = BFGSopt(theta0,H0,false);
      H0Opt = H0;
      theta0Opt = theta0;
      
      
      spc(k,0) = ftrial;
      spc(k,1) = -yopt;
      spc(k,2) = theta0Opt(0);
      spc(k,3) = theta0Opt(1);
      
      ftrial -=fdelta;
      k++;
      if(k>=fstep) break;
    }
    
    return spc;
}



//BFGS core optimization
// pass value by reference
double simple_gpModel::BFGSopt(vec &theta0, mat &H0,
                               bool identity){
  vec deltaFk, deltaFkp1, Xk, Xkp1, pk, sk, yk;
  mat Hk,tmpOut;
  double epsilon, rho, Fk,Fkp1;
  double armijoRight, curvL,curvR, C1, C2, tmpIn;
  double beta1,stepsize,alow,ahigh;
  int iter,iterMax;
  //initial value
  iter = 0;
  iterMax = 1000;
  C1 = 0.01;
  C2 = 0.8;
  beta1 = 0.6;
  epsilon = 1e-4;
  Hk = H0;
  Xk = theta0;
  stepsize = 1;
  // the covariance matrix, only once
  Fk = gp_mLoglik(Xk,true,true);
  if(!flag_inv) return largeV;  
  Fkp1 = Fk;
  deltaFk = gp_DmLoglik(Xk);
  
  while(norm(deltaFk)>epsilon && iter<iterMax){
    iter++;
    // search direction
    pk = -Hk*deltaFk;
    // Wolfe Condition
    stepsize = 1;
    alow = 0;
    ahigh = 10.0;
    tmpIn = dot(deltaFk,pk);
    while(true){
      if(ahigh-alow<1e-3) stepsize = alow;
      Xkp1 = Xk + stepsize*pk;
      Fkp1 = gp_mLoglik(Xkp1,true,true);
      
      //check the armijo condition
      armijoRight =  Fk + C1*stepsize*tmpIn;
      if(Fkp1 > armijoRight || !flag_inv){
        ahigh = stepsize;
        beta1 = 0.66;
      }
      else{
        // check the curvature condition
        deltaFkp1 = gp_DmLoglik(Xkp1);
        if(ahigh-alow<0.0001) break;
        curvR = tmpIn;
        curvL = dot(deltaFkp1,pk);
        if(fabs(curvL)<fabs(C2*curvR)){
          break;
        }else if(curvL > fabs(curvR)){
          ahigh = stepsize;
          beta1 = 0.66;
        }else{
          alow = stepsize;
          beta1 = 0.33;
        }
      }
      if(ahigh-alow<1e-3) break;
      //new try step
      stepsize = alow + (ahigh-alow) * beta1;
    }
    //after fixing alpha
    sk = Xkp1 - Xk;
    if(!flag_inv) return largeV;  
    deltaFkp1 = gp_DmLoglik(Xkp1);
    yk = deltaFkp1 - deltaFk;
    
    // update H
    rho = 1/dot(yk,sk);
    tmpOut = eye<mat>(2,2) - rho*sk*yk.t();
    //	if(identity && iter==1)
    //    Hk = rho/dot(yk,yk) * eye<mat>(2,2);
    Hk = tmpOut *Hk * tmpOut.t() +rho*sk*sk.t();
    // k = k + 1
    Xk = Xkp1;
    Fk = Fkp1;
    deltaFk = deltaFkp1;
  }
  theta0 = Xk;
  H0 = Hk;
  return Fkp1;
}






RCPP_MODULE(varStar_m3){
  class_<varStar>("varStar")
  .constructor<NumericVector,
  NumericVector,
  NumericVector>()
  ;    
  class_ <simple_gpModel>( "simple_gpModel")
    .derives<varStar>("varStar")
    .constructor<NumericVector, NumericVector,
  NumericVector>()
    .method("set_theta", &simple_gpModel::set_theta)
    .method("set_freq", &simple_gpModel::set_freq)
    .method("gp_mLoglik", &simple_gpModel::gp_mLoglik)
    .method("predict", &simple_gpModel::predict)
    .method("BFGSopt", &simple_gpModel::BFGSopt)
    .method("freq_est", &simple_gpModel::freq_est)
  ;
}// MODULE




