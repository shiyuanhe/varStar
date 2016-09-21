
#ifndef MODELGPCLASS2
#define MODELGPCLASS2
#define ARMA_DONT_PRINT_ERRORS

#include "model.h"

class simple_gpModel: public varStar{
protected:
  // prior hyper-para
  
  //gp_m: prior overall mean of the galaxy
  //sigmamSq: prior variance of the stars
  //sigmabSq: prior for beta1 and beta2
  // mu1, mu2: prior mean for theta1 and theta2
  //sigmaSq1, sigmaSq2: prior variance
  double gp_m, sigmamSq, sigmabSq;
  double mu1, mu2,sigmaSq1, sigmaSq2;
  bool flag_inv;
  
  double fmin, fmax, fdelta;
  //gp_y: (=mag-gp_m) de-prior-meaned version of mag, 
  //theta: parameters of kernels
  vec gp_y, theta;
  
  
  //new obs point for prediction
  vec Xnew;
  
  //Ht(H^T): design matrix (1 cos sin),
  //HtBH: uncertainty measurement
  //HtBH = H^T * Sigma_gamma * H
  //Sigma_gamma = diag(sigmamSq, sigmabSq, sigmabSq) 
  //Ht_star: design matrix for new data
  vec gamma0;
  mat Ht, HtBH,Sigma_gamma;
  
  //gp_Ky : all kernel including basis uncertainty and
  //         measurement error
  //        gp_Kc + HtBH
  //its inverse, cholesky factor and determinant
  mat gp_Ky, gp_Ky_inv,gp_Ky_chol;
  double gp_Ky_logdet;
  vec alpha; // gp_Ky_inv * gp_y
  
  //gp_Kc: kernel with only gp and measurement error
  //gp_K21: K(new, train), covariance only involving kernels
  //gp_K22: K(new, new)
  mat gp_K21, gp_K22, gp_Kc, gp_Kc_inv;
  
  
  void compute_CovM(vec &, vec &, int, bool,bool);
  mat gp_get_Ht(vec );
  
  //result of kenrel (and its derivative) computation
  mat kernel_K, kernel_partialK;
  //squared exponential
  void kernel_SE(int, int, int , vec &, vec & ,bool);
  
  
  double largeV;
public:
  
  double gp_mLoglik(vec, bool, bool); //minus log-lik
  vec gp_DmLoglik(vec); // deriv of the above
  mat freq_est();
  double BFGSopt(vec&, mat&, bool);
  
  

  // the log-likelihood for the white noise model
  //double gp_mLoglik_wn();
  
  void set_freq(double, double);
  void set_theta(NumericVector);
  
  List predict(NumericVector);
  
  //MJD, mag, sigma, prior hyperparameters
  simple_gpModel(NumericVector, NumericVector,
                 NumericVector);
};


#endif
