#define ARMA_DONT_PRINT_ERRORS
#include <cstdlib>
#include <string>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>

using namespace std;
using namespace arma;

// These three functions are dealing with file names or read lines
string rm_front_spaces(string in_string) {
  string out_string = in_string;
  while (out_string.substr(0,1) == " ") {
    out_string.erase(0,1);
  }
  return(out_string);
}

string trim_path_from_file_name(string in_string) {
  string out_string = in_string;
  size_t slash_pos = in_string.find_last_of("/\\");
  return(out_string.substr(slash_pos+1));
}

string rm_last_if_slash(string in_string) {
  string out_string = in_string;
  if (in_string.substr(in_string.size()-1) == "/")
    out_string = out_string.substr(0,in_string.size()-1);
  return(in_string);
}

// Define a class for the calculation to avoid many parameter passes
class GP_Model{
protected:
  const double PI = 3.1415926536,
    epsilon = 1e-4,
    C1 = 0.001,
    C2 = 0.8;
  double gp_m, // mean magnitude of the star
    sigmamSq, // variance of the mean magnitude. It will be an empirical constant, say, 2
    sigmabSq; // variance of the periodic component. It will be also an empirical constant, 2
  bool flag_inv; // whether a matrix is invertible or not

  mat spc, spc_1, spc_2; // short for spectrum. Five columns: frequency, log-likelihood, optimal theta1, optimal theta2, curvature of log-likelihood at optimal theta1 and theta2

  vec gp_y, theta; // magnitude - mean magnitude, theta = { theta_1, theta_2)

  vec gamma0;
  mat Ht, // Transpose of the design matrix for the periodic component
    Sigma_gamma, // Covariance matrix ~ diag(sigmamSq, sigmabSq, sigmabSq)
    HtBH; // = Ht * Sigma_gamma * H
  mat gp_get_Ht(vec); // method to get the design matrix

  void compute_CovM(vec &, int); // compute the kernel
  
  mat gp_Kc, // square exp kernel + observation uncertainties
    gp_Ky, // Combined kernel
    gp_Ky_inv, // inverse of the above
    gp_Ky_chol; //Cholesky factorization
  double gp_Ky_logdet; // log (|Ky|).  LOG_e not LOG_10, unfortunately...
  vec alpha; // = gp_Ky_inv * gp_y, for calculation convenience

  mat kernel_K, // Used to store the results of the square exponential kernel
    kernel_partialK; // and its partial derivative (First order)
  void kernel_SE(int, int, int, vec &); // method to compute the above
  
  double largeV; // A very large value.
  
public:
  vec MJD, // Not necessarily MJD, can be JD - 2450000, or what ever time units. By default they are double precision.
    mag, // magnitudes
    mag_sigma; // uncertainties of the magnitudes
  double freq; // variable to store a trial frequency
  int n_obs, n_trials; // number of observations

  double gp_mLoglik(vec); // calculate negative log-likelihood
  vec gp_DmLoglik(vec); // derivative of the above
  double BFGSopt(vec &, mat &); // optimization over theta_1 and theta_2 using the BFGS method

  void set_freq(double); // set the trial frequency
  void freq_est(vec, int); // loop over all the trial frequencies
  vec check_bad(int);
  vec check_breaks(double);
  void load_input(string, string); // input file name, delimiter. (3-column file with mjd, mag, err)
  void spc_output(string); //Write the result to a file
};



// Calculate the square exponential kernel
void GP_Model::kernel_SE(int p1, int p2, int k, vec & t) {
  int i, j;
  double tmpE, tmpE_deriv, tmpK, tmpDiff, eTheta1, eTheta2;
  kernel_K = mat(n_obs, n_obs, fill::zeros);
  kernel_partialK = kernel_K;
  eTheta1 = exp(theta(p1));
  eTheta2 = exp(theta(p2));
  for (i = 0; i < n_obs; i++) {
    for (j = i; j < n_obs; j++) {
      tmpDiff = t(i) - t(j);
      tmpE = tmpDiff * tmpDiff / (2 * eTheta2 * eTheta2);
      if (k == 0) {
	tmpK = eTheta1 * eTheta1 * exp(-tmpE);
	kernel_K(i, j) = tmpK;
	kernel_K(j, i) = tmpK;
      } else if (k == 1) {
	tmpK = 2 * eTheta1 * exp(-tmpE) * eTheta1;
	kernel_partialK(i, j) = tmpK;
	kernel_partialK(j, i) = tmpK;
      } else if (k == 2) {
	tmpK = eTheta1 * eTheta1 * exp(-tmpE);
	tmpE_deriv = 2 * tmpE;
	tmpK *= tmpE_deriv;
	kernel_partialK(i, j) = tmpK;
	kernel_partialK(j, i) = tmpK;
      } else {
	cout << " >> PLEASE CHECK Kernel_SE--k" << endl;
      }
    }
  }
}


// Computer the combined kernel, or covariance matrix
void GP_Model::compute_CovM(vec & t, int k) {
  kernel_SE(0, 1, 0, t);
  gp_Kc = kernel_K;
  for (int i = 0; i < n_obs; i++) {
    gp_Kc(i, i) += mag_sigma(i) * mag_sigma(i) + 0.01;
  }
  gp_Ky = gp_Kc + HtBH; // sigmamSq and sigmabSq are in HtBH
  flag_inv = chol(gp_Ky_chol, gp_Ky);
  if (flag_inv) {
    vec gp_Ky_diag = gp_Ky_chol.diag();
    gp_Ky_logdet = sum(log(gp_Ky_diag));
    gp_Ky_logdet = 2 * gp_Ky_logdet;
    alpha = solve(trimatl(gp_Ky_chol.t()), gp_y);
    alpha = solve(trimatu(gp_Ky_chol), alpha);
  }
}

// Compute Q
double GP_Model::gp_mLoglik(vec theta_) {
  double mLoglik;
  theta = theta_;
  compute_CovM(MJD, 0);
  if (flag_inv) {
    mLoglik = 0.5 * (as_scalar(gp_y.t() * alpha) + gp_Ky_logdet + n_obs * log(2 * PI));
  } else {
    mLoglik = largeV;
  }
  return mLoglik;
}

// Compute partial derivatives of Q
vec GP_Model::gp_DmLoglik(vec theta_) {
  vec partialDeri(2, fill::zeros);
  theta = theta_;
  mat tempI = eye(n_obs, n_obs);
  gp_Ky_inv = solve(trimatl(gp_Ky_chol.t()), tempI);
  gp_Ky_inv = solve(trimatu(gp_Ky_chol), gp_Ky_inv);
  mat tempAATK = alpha * alpha.t() - gp_Ky_inv;
  for (int j = 0; j < 2; j++) {
    kernel_SE(0, 1, j+1, MJD);
    partialDeri(j) = -0.5 * trace(tempAATK * kernel_partialK);
  }
  return partialDeri;
}

// optimization with quasi-Newton method
double GP_Model::BFGSopt(vec & theta0, mat & H0) {
  vec deltaFk, deltaFkp1, Xk, Xkp1, pk, sk, yk;
  mat Hk,tmpOut;
  double rho, Fk,Fkp1;
  double armijoRight, curvL,curvR, C1, C2, tmpIn;
  double beta1,stepsize,alow,ahigh;
  int iter,iterMax;
  iter = 0;
  iterMax = 50;
  C1 = 0.01;
  C2 = 0.8;
  beta1 = 0.6;
  Hk = H0;
  Xk = theta0;
  stepsize = 1;
  // the covariance matrix, only once
  Fk = gp_mLoglik(Xk);
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
      Fkp1 = gp_mLoglik(Xkp1);
	    
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
	}else if(curvL > abs(curvR)){
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

void GP_Model::set_freq(double f) {
  freq = f;
  Ht = gp_get_Ht(MJD);
  HtBH = Ht * Sigma_gamma * Ht.t();
}

mat GP_Model::gp_get_Ht(vec t) {
  mat Ht_tmp = mat(n_obs, 3, fill::ones);
  vec phase = 2 * PI * freq * t;
  Ht_tmp(span::all, 1) = cos(phase);
  Ht_tmp(span::all, 2) = sin(phase);
  return Ht_tmp;
}

void GP_Model::freq_est(vec all_freqs, int restart_sep) {
  n_trials = all_freqs.size();
  double y_min, y_opt;
  bool restart;
  int i, j, k,
    n_trial_theta = 66;
  spc = mat(n_trials, 8, fill::zeros);
  spc_1 = spc;
  spc_2 = spc;
  vec theta0 = vec(2, fill::zeros), theta0_opt, theta0_fix;
  mat H0, H0_opt, H0_pre, H0_fix,
    theta_guess = mat(n_trial_theta, 2, fill::randu);
  theta0_fix = theta0;
  H0_fix = eye<mat> (2, 2);

  theta_guess = (theta_guess - 0.5) * 20;
  theta_guess(0, 0) = 1;
  theta_guess(0, 1) = 1;
  theta_guess(1, 0) = 8;
  theta_guess(1, 1) = 1;
  theta_guess(2, 0) = 1;
  theta_guess(2, 1) = 8;
  theta_guess(3, 0) = 8;
  theta_guess(3, 1) = 8;
  theta_guess(4, 0) = -0.5;
  theta_guess(4, 1) = 4;
  theta_guess(5, 0) = -2.5;
  theta_guess(5, 1) = 5;
    
  for (k = 0; k < n_trials; k++) {
    if ((k % restart_sep) == 0) {
      restart = true;
    } else {
      restart = false;
    }
    set_freq(all_freqs(k));
    if (restart) {
      y_opt = largeV;
      for (i = 0; i < n_trial_theta; i ++) {
	theta0(0) = theta_guess(i, 0);
	theta0(1) = theta_guess(i, 1);
	H0 = eye<mat> (2, 2); // Identity Maxtrix as initial model Hessian
	y_min = BFGSopt(theta0, H0);
	if (y_min < y_opt) {
	  y_opt = y_min;
	  H0_opt = H0;
	  theta0_opt = theta0;
	}
      }
    }
    y_min = BFGSopt(theta0_opt, H0_opt);
    spc_1(k, 0) = freq;
    spc_1(k, 1) = -1 * y_min;
    spc_1(k, 2) = theta0_opt(0);
    spc_1(k, 3) = theta0_opt(1);
    spc_1(k, 4) = H0_opt(0, 0);
    spc_1(k, 5) = H0_opt(0, 1);
    spc_1(k, 6) = H0_opt(1, 0);
    spc_1(k, 7) = H0_opt(1, 1);
  }
  for (k = n_trials - 1; k >= 0; k--) {
    if ((k % restart_sep) == 0) {
      restart = true;
    } else {
      restart = false;
    }
    set_freq(all_freqs(k));
    if (restart) {
      theta0_opt(0) = spc_1(k, 2);
      theta0_opt(1) = spc_1(k, 3);
      H0_opt(0, 0) = spc_1(k, 4);
      H0_opt(0, 1) = spc_1(k, 5);
      H0_opt(1, 0) = spc_1(k, 6);
      H0_opt(1, 1) = spc_1(k, 7);
    } else if (k == n_trials - 1) {
      y_opt = largeV;
      for (i = 0; i < n_trial_theta; i ++) {
	theta0(0) = theta_guess(i, 0);
	theta0(1) = theta_guess(i, 1);
	H0 = eye<mat> (2, 2);
	y_min = BFGSopt(theta0, H0);
	if (y_min < y_opt) {
	  y_opt = y_min;
	  H0_opt = H0;
	  theta0_opt = theta0;
	}
      }
    }
    y_min = BFGSopt(theta0_opt, H0_opt);
    spc_2(k, 0) = freq;
    spc_2(k, 1) = -1 * y_min;
    spc_2(k, 2) = theta0_opt(0);
    spc_2(k, 3) = theta0_opt(1);
    spc_2(k, 4) = H0_opt(0, 0);
    spc_2(k, 5) = H0_opt(0, 1);
    spc_2(k, 6) = H0_opt(1, 0);
    spc_2(k, 7) = H0_opt(1, 1);
  }
  
  for (k = 0; k < n_trials; k++) {
    if (spc_1(k, 1) >= spc_2(k, 1)) {
      spc(k, 0) = spc_1(k, 0);
      spc(k, 1) = spc_1(k, 1);
      spc(k, 2) = spc_1(k, 2);
      spc(k, 3) = spc_1(k, 3);
      spc(k, 4) = spc_1(k, 4);
      spc(k, 5) = spc_1(k, 5);
      spc(k, 6) = spc_1(k, 6);
      spc(k, 7) = spc_1(k, 7);
    } else {
      spc(k, 0) = spc_2(k, 0);
      spc(k, 1) = spc_2(k, 1);
      spc(k, 2) = spc_2(k, 2);
      spc(k, 3) = spc_2(k, 3);
      spc(k, 4) = spc_2(k, 4);
      spc(k, 5) = spc_2(k, 5);
      spc(k, 6) = spc_2(k, 6);
      spc(k, 7) = spc_2(k, 7);
    }
  }

  // Check bad points and fix the initial guess of H0
  int last_idx_size = -1;
  for (j = 0; j < 30; j++) {
    double diff_val = 0.2;
    int inc_length = 1;
    vec idx = check_breaks(diff_val);
    vec idx_2 = check_breaks(diff_val*10);
    vec idx_refit = vec(idx.size() * (2*inc_length + 1));
    int jk = 0, ji, jj, max_idx = 0;

    if (last_idx_size <= int (idx.size())) {
      break;
    }
    last_idx_size = idx.size();

    if (idx_2.size() > 10)
      break;

    if (idx.size() > 5 && max(spc(span(0, n_trials-1),1)) > 0) 
      break;

    if (idx.size() > 0) {
      for (int ji = 0; ji < int (idx.size()); ji++) {
	for (jj = -inc_length; jj <= inc_length; jj++) {
	  idx_refit(jk) = idx(ji) + jj;
	  jk++;
	}
      }
      idx_refit = unique(idx_refit);
      int n_refit = idx_refit.size();
      cout << "  Refit for the " << j+1 << "th iteration:";
      cout << "   Fixing at " << n_refit << " positions..." << endl;

      H0_fix(0, 0) = median(spc(span(0, n_trials-1),4));
      H0_fix(0, 1) = median(spc(span(0, n_trials-1),5));
      H0_fix(1, 0) = median(spc(span(0, n_trials-1),6));
      H0_fix(1, 1) = median(spc(span(0, n_trials-1),7));
      theta_guess(0, 0) = median(spc(span(0, n_trials-1),2));
      theta_guess(0, 1) = median(spc(span(0, n_trials-1),3));
      theta_guess(1, 0) = (randu()-0.5) * j;
      theta_guess(1, 1) = (randu()-0.5) * j;
      
      vec tmp_vec = spc(span::all,1);
      double max_val = max(tmp_vec);
      for (int i_max = 0; i_max < n_trials; i_max++) {
	if (tmp_vec(i_max) == max_val) 
	  max_idx = i_max;
      }
      theta_guess(2, 0) = spc(max_idx, 2);
      theta_guess(2, 1) = spc(max_idx, 3);
      
    
      for (jk = 0; jk < n_refit; jk++) {
	if (idx_refit(jk) >= 0 && idx_refit(jk) < n_trials) {
	  set_freq(all_freqs(idx_refit(jk)));
	  if (jk >= 0) {
	    y_opt = largeV;
	    y_min = largeV;
	    for (ji = 0; ji < 3; ji++) {
	      theta0(0) = theta_guess(ji, 0);
	      theta0(1) = theta_guess(ji, 1);
	      H0 = H0_fix; 
	      y_min = BFGSopt(theta0, H0);
	      if (y_min < y_opt) {
		y_opt = y_min;
		H0_opt = H0; // Find a good Hessian Matrix
		theta0_opt = theta0;
	      }
	    }
	    if (idx_refit(jk) > 1) {
	      theta0(0) = spc(idx_refit(jk)-1,2);
	      theta0(1) = spc(idx_refit(jk)-1,3);
	    }
	    y_min = BFGSopt(theta0, H0_opt);
	    if (y_min < y_opt) {
	      y_opt = y_min;
	      H0_opt = H0;
	      theta0_opt = theta0;
	    }
	  }
	  y_min = BFGSopt(theta0_opt, H0_opt);
	  spc(idx_refit(jk), 0) = freq;
	  spc(idx_refit(jk), 1) = -1 * y_min;
	  spc(idx_refit(jk), 2) = theta0_opt(0);
	  spc(idx_refit(jk), 3) = theta0_opt(1);
	  spc(idx_refit(jk), 4) = H0_opt(0, 0);
	  spc(idx_refit(jk), 5) = H0_opt(0, 1);
	  spc(idx_refit(jk), 6) = H0_opt(1, 0);
	  spc(idx_refit(jk), 7) = H0_opt(1, 1);
	} 
      }
    } else {
      break;
    }
  }
}

vec GP_Model::check_breaks(double diff_val) {
  vec tmp_idx(n_trials), tmp_diff;
  int k, j = 0;
  tmp_diff = spc(span(2, n_trials-1), 1) - 2 * spc(span(1, n_trials-2), 1) + spc(span(0, n_trials-3), 1);
  tmp_diff = abs(tmp_diff);
  for (k = 0; k < n_trials - 2; k++) {
    if (tmp_diff(k) > diff_val) {
      tmp_idx(j) = k;
      j++;
    }
  }
  j--;
  if (j > 0) {
    tmp_idx = tmp_idx.head(j);
    return tmp_idx;
  } else {
    vec no_idx;
    return no_idx;
  }
}

vec GP_Model::check_bad(int col_idx) {
  vec tmp_idx(n_trials);
  double tmp_mean, tmp_sd;
  int j=0;
  tmp_mean = mean(spc(span(0,n_trials-1),col_idx));
  tmp_sd = stddev(spc(span(0,n_trials-1),col_idx));
  for (int k = 0; k < n_trials; k++) {
    if (spc(k, col_idx) > (tmp_mean + 2.3*tmp_sd) || spc(k, col_idx) < (tmp_mean - 2.3*tmp_sd)) {
      tmp_idx(j) = k;
      j++;
    }
  }
  j--;
  tmp_idx = tmp_idx.head(j);
  return tmp_idx;
}

void GP_Model::load_input(string f_lc, string delimiter) {
  string line, firstline, tsstring, rest_string;
  int header, isub, iline;
  // [[1]] load light curve file (1): Count the number of observations. Autodetect file headers
  header = 0;
  ifstream cntmyfile(f_lc);
  if (cntmyfile.is_open()) {
    getline(cntmyfile,firstline);
    rest_string = rm_front_spaces(firstline);
    header = 1;
    if (rest_string.substr(0,1) == "-")
      header = 0;
    for (isub=0;isub<=9;isub++) {
      tsstring = to_string(isub);
      if (rest_string.substr(0,1) == tsstring)
	header = 0;
    }
    cntmyfile.unsetf(ios_base::skipws);
    n_obs = count(istream_iterator<char>(cntmyfile), istream_iterator<char>(), '\n');
    if (header == 0)
      n_obs = n_obs + 1;
    cntmyfile.close();
  } else {
    cout << "Cannot open file " << f_lc << "\n";
    abort();
  }
  // [[2]] Load light curves (2): Read file into vectors MJD, mag , mag_sigma, weight
  MJD = vec(n_obs);
  mag = vec(n_obs);
  mag_sigma = vec(n_obs);
  ifstream myfile(f_lc);
  if (myfile.is_open()) {
    if (header == 1)
      getline(myfile,line);
    for (iline=0; iline < n_obs; iline++) {
      getline(myfile,line);
      rest_string = rm_front_spaces(line);
      MJD[iline] = stod(rest_string);
      rest_string = rm_front_spaces(rest_string.substr(rest_string.find(delimiter),rest_string.length()));
      mag[iline] = stod(rest_string);
      rest_string = rm_front_spaces(rest_string.substr(rest_string.find(delimiter),rest_string.length()));
      mag_sigma[iline] = stod(rest_string);
    }
    myfile.close();
  } else {
    cout << "Cannot open file " << f_lc << "\n";
    abort();
  }

  //[[2.5]] A simple 3-sigma clip to remove significant outliers
  float mag_stddev = stddev(mag);
  float mag_mean = mean(mag);
  float outlier_sigma = 3;
  int kpt_n_obs = 0, kpt_index = 0;
  vec kpt_mag, kpt_MJD, kpt_mag_sigma;
  for (int i_foo=0; i_foo<n_obs; i_foo++) {
    if (abs(mag[i_foo] - mag_mean) < outlier_sigma*mag_stddev) {
      kpt_n_obs ++;
    }
  }
  kpt_mag = vec(kpt_n_obs);
  kpt_MJD = vec(kpt_n_obs);
  kpt_mag_sigma = vec(kpt_n_obs);
  for (int i_foo=0; i_foo<n_obs; i_foo++) {
    if (abs(mag[i_foo] - mag_mean) < outlier_sigma*mag_stddev) {
      kpt_mag[kpt_index] = mag[i_foo];
      kpt_MJD[kpt_index] = MJD[i_foo];
      kpt_mag_sigma[kpt_index] = mag_sigma[i_foo];
      kpt_index ++;
    }
  }
  mag = kpt_mag;
  MJD = kpt_MJD;
  mag_sigma = kpt_mag_sigma;
  n_obs = kpt_n_obs;

  // [[3]] Initialize some parameters
  gp_m = mean(mag);
  gp_y = mag - gp_m;
  gamma0 = vec(3,fill::zeros);
  gamma0(0) = gp_m;
  sigmamSq = 5;
  sigmabSq = 5;
  Sigma_gamma = mat(3, 3, fill::zeros);
  Sigma_gamma(0,0) = sigmamSq;
  Sigma_gamma(1,1) = sigmabSq;
  Sigma_gamma(2,2) = sigmabSq;
  largeV = 1.0e30;
}

void GP_Model::spc_output(string f_output) {
  ofstream myfile (f_output);
  if (myfile.is_open()){
    for (int i = 0; i < n_trials; i++)
      {
	if (spc(i,0) != 0) {
	  myfile << fixed << setprecision(6) << spc(i,0)
		 << "  " << spc(i,1)
		 << "  " << spc(i,2)
		 << "  " << spc(i,3)
		 // << "  " << spc(i,4)
		 // << "  " << spc(i,5)
		 // << "  " << spc(i,6)
		 // << "  " << spc(i,7)
		 <<endl;
	}
      }
    myfile.close();
  } 
}

int main(int argc, char* argv[ ]) {


  string delimiter = " "; // delimiter for the input light curve file. " " is good enough for multiple-space delimited files
  string output_dir = "./gp_spectra/"; // directory that used to save output files (frequency spectra)
  string output_extension = ".gp.dat"; // extension of output file names.
  string model_name = argv[0];
  string syntax_prom = "********** Syntax ***********\n " 
    + model_name + 
    " -f light_curve_file_name (Compute a single light curve) \n OR \n "
    + model_name + 
    " -l list_file_name (Compute a list of light curves)\n*****************************\n";

  int restart_sep = 50;
  vec all_freqs, freqs_hp, freqs_lp;
  double p1 = 2000.;
  double p2 = 100.;
  double p3 = 10.;
  double hp_resolution = 1.e-5; // in frequency
  double lp_resolution = 0.2; // in days
  int n_freqs_hp, n_freqs_lp;
  n_freqs_hp = int((1./p2 - 1./p1) / hp_resolution);
  n_freqs_lp = int((p2 - p3) / lp_resolution);
  
  all_freqs = vec(n_freqs_hp + n_freqs_lp);
  
  for (int i = 0; i < n_freqs_hp; i++) {
    all_freqs(i) = 1./p2 - double(i) * hp_resolution;
  }

  for (int i = 0; i < n_freqs_lp; i++) {
    all_freqs(i + n_freqs_hp) = 1./(p2 - double(i) * lp_resolution);
  }
  all_freqs = unique(all_freqs);
  
  
  if (argc > 2) {
    string f_lc, f_list, file_type = argv[1];
    string f_output, f_lc_trim, shell_command;
      
    output_dir = rm_last_if_slash(output_dir);
    shell_command = "mkdir -p " + output_dir;
    const char *ts_command = shell_command.c_str();
    system(ts_command);
    if (file_type == "-f") {
      f_lc = argv[2];
      GP_Model lc;
      lc.load_input(f_lc, delimiter);		       
      lc.freq_est(all_freqs, restart_sep);
      f_lc_trim = trim_path_from_file_name(f_lc);
      f_output = output_dir + "/" + f_lc_trim + output_extension;
      lc.spc_output(f_output);
    } else if (file_type == "-l") {
      f_list = argv[2];
      ifstream list_file (f_list);
      if (list_file.is_open()) {
	unsigned gp_counter = 0;
	GP_Model lc;
	while (!list_file.eof()) {
	  list_file >> f_lc;
	  lc.load_input(f_lc, delimiter);
	  lc.freq_est(all_freqs, restart_sep);
	  f_lc_trim = trim_path_from_file_name(f_lc);
	  f_output = output_dir + "/" + f_lc_trim + output_extension;
	  lc.spc_output(f_output);
	  gp_counter++;
	  string f_status = "status/" + f_list + ".status.txt";
	  f_status.replace(f_status.begin()+7,f_status.begin()+12,"");
	  ofstream out_status(f_status);
	  out_status << gp_counter << "   " << f_lc_trim << "\n";
	  out_status.close();
	  cout << gp_counter << "   " << f_lc_trim << "\n";
	}		
	list_file.close();
      }
      else {
	cout << "Cannot open file " << f_list << ". Please double check the file name.\n";
	abort();
      }
    } else {
      cout << syntax_prom;
    } 
  } else {
    cout << syntax_prom;
  } 
  return 0; 
}
// v = 2.8
