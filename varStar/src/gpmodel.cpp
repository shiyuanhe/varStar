#include "gpmodel.h"

void gpModel::kernel_SE(int p1, int p2, int k){
    int i,j;
    double tmp, tmpK, tmpDiff;
    kernel_K = mat(nK1, nK2, fill::zeros);
    for(i=0; i<nK1; i++){
        for(j=0; j<nK2; j++){
            tmpDiff = gp_y1(i) - gp_y2(j);
            tmpDiff = tmpDiff * tmpDiff;
            // if(i==0 && j==1) Rcout<<tmpDiff<<std::endl;
            //long term trend
            tmp = tmpDiff/(2*theta(p2)*theta(p2));
            //if(i==0 && j==1) Rcout<<tmp<<std::endl;
            if(k==0){
                tmpK = theta(p1)*theta(p1)*exp(-tmp);
                //if(i==0 && j==1) Rcout<<tmpK<<std::endl;
                kernel_K(i,j) = tmpK;
            }else if(k==1){
                tmpK = 2*theta(p1)*exp(-tmp);
                kernel_partialK(i,j) = tmpK;
            }else if(k==2){
                tmpK = theta(p1)*theta(p1)*exp(-tmp);
                tmp = 2*tmp/theta(p2);
                tmpK *= tmp;
                kernel_partialK(i,j) = tmpK;
            }
        }
    }
}

void gpModel::kernel_periodic(int p1, int p2, int p3, int k){
    int i, j;
    double tmp, tmpK, tmpDiff, sinDiff;
    kernel_K = mat(nK1, nK2, fill::zeros);
    for(i=0; i<nK1; i++){
        for(j=0; j<nK2; j++){
            tmpDiff = gp_y1(i) - gp_y2(j);
            
            sinDiff = sin(PI * tmpDiff * freq);
            tmpDiff = tmpDiff * tmpDiff;
            sinDiff = sinDiff * sinDiff;
            //long term trend
            tmp = tmpDiff/(2*theta(p2)*theta(p2));
            tmp += 2*sinDiff/(theta(p3)*theta(p3));
            //if(i==0 && j==1) Rcout<<tmp<<std::endl;
            if(k==0){
                tmpK = theta(p1)*theta(p1)*exp(-tmp);
                //if(i==0 && j==1) Rcout<<tmpK<<std::endl;
                kernel_K(i,j) = tmpK;
            }else if(k==1){
                tmpK = 2*theta(p1)*exp(-tmp);
                kernel_partialK(i,j)= tmpK;
            }else if(k==2){
                tmpK = theta(p1)*theta(p1)*exp(-tmp);
                tmp = tmpDiff/(theta(p2)*theta(p2)*theta(p2));
                tmpK *= tmp;
                kernel_partialK(i,j) = tmpK;
            }else if(k==3){
                tmpK = theta(p1)*theta(p1)*exp(-tmp);
                tmp = 4*sinDiff/(theta(p3)*theta(p3)*theta(p3));
                tmpK *= tmp;
                kernel_partialK(i,j) = tmpK;
            }
        }
    }
}


//compute different part of variance/covariance matrix
//task==0 and nK1==nK2, is to compute training data
//assign values to gp_y1, gp_y2, nK1, nK2 before computing
//task==0: compute all
void gpModel::compute_CovM(int taskI){
    gp_K = mat(nK1, nK2, fill::zeros);
    
    if(taskI==0 || taskI==1){
        kernel_SE(0, 1, 0);
        gp_K = gp_K + kernel_K;
        
    }
    
    if(taskI==0 || taskI==2){
        kernel_periodic(2,3,4,0);
        gp_K = gp_K + kernel_K;
        
    }
    
    if(taskI==0 || taskI==3){
        kernel_SE(5, 6, 0);
        gp_K = gp_K + kernel_K;
        
    }
    
    if(taskI==0 && nK1==nK2){
        for(int i=0; i<nObs; i++){
            //if(i < 5) Rcout<<1/(weight(i)*weight(i))<<std::endl;
            gp_K(i, i)  += mag_sigma(i)*mag_sigma(i);  //relaxed
        }
        gp_Kinv = gp_K.i();
    }
}


//copmute minus log likelihood of the data
double gpModel::gp_mLoglik(NumericVector theta_){
    double mLoglik = 0, signL;
    theta = vec(theta_);
    compute_CovM(0);
    
    log_det(mLoglik, signL, gp_K);
    
    mLoglik += as_scalar(gp_y.t()*gp_Kinv*gp_y);
    return mLoglik;
}

//????????? the minus sign
//compute the partial derivative of the minus log likelihood
vec gpModel::gp_DmLoglik(NumericVector theta_){
    vec partialDeri(7);
    theta = vec(theta_);
    compute_CovM(0);
    colvec alpha;
    mat prodTmp;
    alpha = gp_Kinv * gp_y;
    prodTmp = alpha*alpha.t() - gp_Kinv;
    //long term
    kernel_SE(0,1,1);
    partialDeri(0) = trace(prodTmp*kernel_partialK);
    kernel_SE(0,1,2);
    partialDeri(1) = trace(prodTmp*kernel_partialK);
    
    //periodic
    kernel_periodic(2,3,4,1);
    partialDeri(2) = trace(prodTmp*kernel_partialK);
    kernel_periodic(2,3,4,2);
    partialDeri(3) = trace(prodTmp*kernel_partialK);
    kernel_periodic(2,3,4,3);
    partialDeri(4) = trace(prodTmp*kernel_partialK);
    
    //noise
    kernel_SE(5,6,1);
    partialDeri(5) = trace(prodTmp*kernel_partialK);
    kernel_SE(5,6,2);
    partialDeri(6) = trace(prodTmp*kernel_partialK);
    
    return partialDeri;
}

//set theta before further MLE
void gpModel::gp_setTheta(NumericVector theta_){
    theta = vec(theta_);
    gp_y = mag - gp_m;
    gp_y1 = MJD;
    gp_y2 = MJD;
    nK1 = nObs;
    nK2 = nObs;
    kernel_partialK = mat(nK1, nK2, fill::zeros);
}

List gpModel::gp_predict(NumericVector Xnew_, int comp){
    Xnew = vec(Xnew_);
    
    //covariance and variance of Xnew
    gp_y1 = Xnew;
    nK1 = Xnew.n_elem;
    gp_y2 = Xnew;
    nK2 = Xnew.n_elem;
    //compute_CovM(0);
    //gp_K22 = gp_K;
    //this computation is not correct
    //No Diag is needed, should not compute the inverse
    
    //variance of Xnew
    gp_y2 = MJD;
    nK2 = nObs;
    compute_CovM(comp);
    gp_K21 = gp_K;
    
    //variance of train
    gp_y = mag - gp_m;
    gp_y1 = MJD;
    nK1 = nObs;
    compute_CovM(0); 
    
    //mat covY = gp_K22 - gp_K21*gp_Kinv*gp_K21.t();
    vec covD = vec(1, fill::zeros);//covY.diag();
    vec predY =gp_m + gp_K21 * gp_Kinv * gp_y;
    return List::create(Named("predy")=predY,
                        Named("var")=covD);
}

gpModel::gpModel(NumericVector MJD_,
                 NumericVector mag_, NumericVector error_, double gpM_):
    varStar(MJD_, mag_, error_)
{
    gp_m = gpM_;
}


void gpModel::set_freq(double fff, double shift){
    freq = fff;
    //new design matrix after changing period
    vec tmpMJD;
    tmpMJD = MJD - min(MJD);
    tmpMJD = tmpMJD * freq;
    tmpMJD = floor(tmpMJD);
    nCyc = max(tmpMJD) + 1;
}


List gpModel::get_fake(NumericVector fMJD_, double noiseFactor){
    RNGScope scope;
    vec fMJD(fMJD_),fMJDStraight;
    int nfObs = fMJD.n_elem; // number of fake observation
    int i;
    mat fLcurve(nfObs, 3, fill::zeros);
    
    //make shift
    //calculate phase and cycles
    double deltaMJD = 0;
    if(nCyc>3)
        deltaMJD =  (nCyc-1) / freq;
    else
        deltaMJD = nCyc / freq;
    
    double shiftMax;
    shiftMax = max(MJD) - min(MJD) - (max(fMJD) - min(fMJD));
    if(shiftMax < 1/ freq) shiftMax = 1 / freq;
    
    double shift = R::runif(0,1) / freq;
    //Rcout<<shift<<std::endl;
    
    vec reflect;
    fLcurve(span::all, 0) = fMJD;
    fMJD = fMJD - min(fMJD) + shift;
    fMJDStraight = fMJD;
    fMJD = fMJD / deltaMJD;
    reflect = floor(fMJD);
    fMJD = fMJD -floor(fMJD);
    int cycP_;
    for(i=0;i<nfObs;i++){
        cycP_ = reflect(i);
        if(cycP_ % 2 ==1) fMJD(i) = 1 - fMJD(i);
    }
    
    fMJD = fMJD*deltaMJD + min(MJD);
    fMJDStraight += min(MJD);
    NumericVector fMJD1 = Rcpp::as<Rcpp::NumericVector>(wrap(fMJD));
    NumericVector fMJD2 = Rcpp::as<Rcpp::NumericVector>(wrap(fMJDStraight));
    List res2 = gp_predict(fMJD2,2);
    List res1 = gp_predict(fMJD1,1);
    List res3 = gp_predict(fMJD1,3);
    
    fLcurve(span::all,1) = vec(as<NumericVector>(res1["predy"]));
    fLcurve(span::all,1) += vec(as<NumericVector>(res3["predy"]));
    fLcurve(span::all,1) += vec(as<NumericVector>(res2["predy"]));
    fLcurve(span::all,1) -= 2*gp_m;
    //add noise level according to mag
    fLcurve(span::all,2) = 0.00078* fLcurve(span::all,1) 
        % fLcurve(span::all,1)-0.02258*fLcurve(span::all,1) 
        + 0.16908;
        fLcurve(span::all, 2) = fLcurve(span::all, 2) * noiseFactor;
        double noiseAdd;
        for(i=0; i<nfObs; i++){
            if(fLcurve(i,1)<15) fLcurve(i,2) = 0.006 * noiseFactor;
            noiseAdd = R::rnorm(0, fLcurve(i,2));
            fLcurve(i,1) += noiseAdd; 
        }
        return List::create(Named("fLcurve") = fLcurve,
                            Named("shift") = shift);
}

RCPP_MODULE(varStar_m2){
    class_<varStar>("varStar")
    .constructor<NumericVector, NumericVector, NumericVector>()
    ;
    
    class_ <gpModel>( "gpModel")
        .derives<varStar>("varStar")
        .constructor<NumericVector, NumericVector, NumericVector, double>()
        .method("set_freq", &gpModel::set_freq)
        .method("gp_setTheta", &gpModel::gp_setTheta)
        .method("gp_mLoglik", &gpModel::gp_mLoglik)
        .method("gp_DmLoglik", &gpModel::gp_DmLoglik)
        .method("gp_predict", &gpModel::gp_predict)
        .method("get_fake", &gpModel::get_fake)
    ;
}
