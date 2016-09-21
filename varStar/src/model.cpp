#include "model.h"
#include <RcppArmadilloExtensions/sample.h>

//In this file, magnitude is modeled as
//a_k*f(phase)+m_k
//f is a function of phase, defined on (0,1)
//a_k and m_k varies from cycle to cycle
//currently impose a_k=1 for all cycles

varStar::varStar(NumericVector MJDI, NumericVector magI,
		    NumericVector errorI){
    MJD = vec(MJDI);
    mag = vec(magI);
    mag_sigma = vec(errorI);
    weight = vec(errorI);
    weight = 1/weight;
    nObs = MJD.n_elem;
}

// Set the frequency
// shift is UNIFORM(0, period)
void varStar::set_freq(double fff, double shift){

}








