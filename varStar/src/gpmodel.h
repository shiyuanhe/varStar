
#ifndef MODELGPCLASS
#define MODELGPCLASS

#include "model.h"

class gpModel: public varStar{
protected:
    int nCyc;
    double gp_m;
    vec gp_y, theta;
    int nK1, nK2;
    vec Xnew, gp_y1, gp_y2;
    mat gp_K, gp_Kinv, gp_K21, gp_K22;
    void compute_CovM(int);
    
    mat kernel_K, kernel_partialK;
    void kernel_SE(int, int, int );//squared exponential
    void kernel_periodic(int p1, int p2, int p3, int k);
    
    
public:
    gpModel(NumericVector, NumericVector, NumericVector, double);
    void set_freq(double, double);
    void gp_setTheta(NumericVector );
    double gp_mLoglik(NumericVector);
    vec gp_DmLoglik(NumericVector);
    List gp_predict(NumericVector, int comp);
    List get_fake(NumericVector fMJD_, double noiseFactor);
};



#endif
