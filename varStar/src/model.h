#ifndef MODELBASECLASS
#define MODELBASECLASS

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <vector>


using namespace Rcpp;
using namespace arma;


class varStar{
protected:
     //input
    //weight is the inverse of input error
    //       is also the inverse of mag_sigma
    vec MJD, mag, weight, mag_sigma;
    double freq;
    int nObs;
        
 public:
    varStar(NumericVector, NumericVector, NumericVector);
    virtual void set_freq(double fff, double shift) ;
    virtual ~varStar() {}
};


#endif //MODELBASECLASS
